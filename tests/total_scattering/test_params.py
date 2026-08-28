"""Tests for locating the facility/instrument configuration.

The reduction historically read its facility configuration straight off the
SNS analysis cluster (``/SNS/NOM/shared/autoreduce/configs``).  That makes the
reduction unrunnable anywhere else - inside a container, on a laptop, or in a
Galaxy job - so the loader has to be able to find the configuration somewhere
else and fall back to a packaged copy.
"""

import json
import os
import unittest
from unittest import mock

from total_scattering import params


class TestParamsLoaderResolution(unittest.TestCase):
    """The config file is looked up through an explicit chain of locations."""

    def setUp(self):
        # Keep every test isolated from whatever the developer has exported.
        self.env_patch = mock.patch.dict(
            os.environ,
            {k: v for k, v in os.environ.items()
             if k not in ("MTS_AUTO_CONFIG_FILE", "MTS_CONFIG_DIR")},
            clear=True,
        )
        self.env_patch.start()
        self.addCleanup(self.env_patch.stop)

    def _write_config(self, directory, **overrides):
        os.makedirs(directory, exist_ok=True)
        contents = {"QParamsProcessing": "0.02,0.002,35.0", "TMIN": 111.0}
        contents.update(overrides)
        path = os.path.join(directory, "auto_config.json")
        with open(path, "w") as handle:
            json.dump(contents, handle)
        return path

    def test_explicit_config_file_wins(self):
        """An explicitly passed file beats every other candidate."""
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            explicit = self._write_config(
                os.path.join(tmp, "explicit"), TMIN=1.0)
            os.environ["MTS_AUTO_CONFIG_FILE"] = self._write_config(
                os.path.join(tmp, "env"), TMIN=2.0)

            loader = params.ParamsLoader("SNS", "NOM", config_file=explicit)

            self.assertEqual(loader.config_file, explicit)
            self.assertEqual(loader.config_params["TMIN"], 1.0)

    def test_env_config_file_used_when_no_explicit_file(self):
        """``MTS_AUTO_CONFIG_FILE`` names a config file directly."""
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            env_file = self._write_config(os.path.join(tmp, "env"), TMIN=2.0)
            os.environ["MTS_AUTO_CONFIG_FILE"] = env_file

            loader = params.ParamsLoader("SNS", "NOM")

            self.assertEqual(loader.config_file, env_file)
            self.assertEqual(loader.config_params["TMIN"], 2.0)

    def test_env_config_dir_used_when_no_config_file(self):
        """``MTS_CONFIG_DIR`` holds an auto_config.json."""
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            config_dir = os.path.join(tmp, "configs")
            self._write_config(config_dir, TMIN=3.0)
            os.environ["MTS_CONFIG_DIR"] = config_dir

            loader = params.ParamsLoader("SNS", "NOM")

            self.assertEqual(loader.config_dir, config_dir)
            self.assertEqual(loader.config_params["TMIN"], 3.0)

    def test_facility_path_used_when_present(self):
        """The on-cluster location beats the packaged copy."""
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            cluster = self._write_config(os.path.join(tmp, "sns"), TMIN=4.0)
            config_loc = {"SNS": {"NOM": cluster}}

            loader = params.ParamsLoader("SNS", "NOM", config_loc_in=config_loc)

            self.assertEqual(loader.config_file, cluster)
            self.assertEqual(loader.config_params["TMIN"], 4.0)

    def test_falls_back_to_packaged_config_off_cluster(self):
        """Off the cluster we still get a usable configuration."""
        config_loc = {"SNS": {"NOM": "/SNS/NOM/does/not/exist/auto_config.json"}}

        loader = params.ParamsLoader("SNS", "NOM", config_loc_in=config_loc)

        self.assertTrue(os.path.isfile(loader.config_file))
        self.assertEqual(
            os.path.dirname(loader.config_file), loader.config_dir)
        # The keys the reduction reads without a default must all be present.
        for key in ("QParamsProcessing", "TMIN", "TMAX", "CacheDir"):
            self.assertIn(key, loader.config_params)

    def test_unknown_instrument_falls_back_to_facility_default(self):
        """An instrument with no packaged config still resolves."""
        config_loc = {"SNS": {"NOSUCH": "/SNS/NOSUCH/nope/auto_config.json"}}

        loader = params.ParamsLoader("SNS", "NOSUCH", config_loc_in=config_loc)

        self.assertTrue(os.path.isfile(loader.config_file))
        self.assertIn("QParamsProcessing", loader.config_params)

    def test_raises_for_unknown_facility(self):
        """An unknown facility is an error, not a silent default."""
        with self.assertRaises(FileNotFoundError):
            params.ParamsLoader("NOWHERE", "NOM", config_loc_in={})


class TestParamsLoaderPackagedDefaults(unittest.TestCase):
    """The packaged configuration has to be shipped, not just importable."""

    def test_packaged_config_dir_exists(self):
        self.assertTrue(os.path.isdir(params.PACKAGED_CONFIG_DIR))

    def test_packaged_nomad_config_is_valid_json(self):
        path = os.path.join(
            params.PACKAGED_CONFIG_DIR, "SNS", "NOM", "auto_config.json")
        self.assertTrue(os.path.isfile(path))
        with open(path) as handle:
            contents = json.load(handle)
        self.assertIsInstance(contents["QParamsProcessing"], str)


if __name__ == '__main__':
    unittest.main()  # pragma: no cover
