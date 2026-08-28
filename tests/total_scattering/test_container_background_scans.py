"""Tests for resolving the container background's scan names.

Every other section of the input JSON accepts either ``Runs`` (run numbers to
be turned into NeXus basenames) or ``Filenames`` (explicit paths). The
container background - ``Sample.Background.Background`` - has to behave the
same way, otherwise a file-based reduction cannot subtract the empty
instrument.
"""

import unittest

import total_scattering.reduction.total_scattering_reduction as ts


class TestContainerBackgroundScans(unittest.TestCase):

    def test_run_numbers_become_nexus_basenames(self):
        sample_background = {"Background": {"Runs": "91607-91608"}}

        scans = ts.container_background_scans(
            sample_background, "NOM", "%s_%d")

        self.assertEqual(scans, "NOM_91607,NOM_91608")

    def test_runs_are_normalised_in_place(self):
        """The reduction reports the expanded run list back to the caller."""
        sample_background = {"Background": {"Runs": "91607-91609"}}

        ts.container_background_scans(sample_background, "NOM", "%s_%d")

        self.assertEqual(
            sample_background["Background"]["Runs"], [91607, 91608, 91609])

    def test_filenames_only_defers_to_the_filenames_override(self):
        """No 'Runs' key is not an error; the file paths are applied later."""
        sample_background = {
            "Background": {"Filenames": ["/data/empty_instrument.nxs"]}}

        scans = ts.container_background_scans(
            sample_background, "NOM", "%s_%d")

        self.assertIsNone(scans)

    def test_no_container_background_section(self):
        self.assertIsNone(
            ts.container_background_scans({}, "NOM", "%s_%d"))

    def test_empty_run_list_is_treated_as_absent(self):
        sample_background = {"Background": {"Runs": []}}

        self.assertIsNone(
            ts.container_background_scans(sample_background, "NOM", "%s_%d"))

    def test_isis_style_file_format(self):
        sample_background = {"Background": {"Runs": "97947"}}

        scans = ts.container_background_scans(
            sample_background, "POLARIS", "%s%d")

        self.assertEqual(scans, "POLARIS97947")


if __name__ == '__main__':
    unittest.main()  # pragma: no cover
