"""Tests for reading the sub-group index that accompanies a facility config.

``group_index.txt`` lives next to ``auto_config.json`` on the analysis cluster
and is only needed when absorption correction regroups detectors.  Off the
cluster it is routinely absent, and that must not abort the reduction.
"""

import os
import tempfile
import unittest

import total_scattering.reduction.total_scattering_reduction as ts


class TestLoadGroupIndex(unittest.TestCase):

    def _write(self, contents):
        handle = tempfile.NamedTemporaryFile(
            "w", suffix=".txt", delete=False)
        handle.write(contents)
        handle.close()
        self.addCleanup(os.remove, handle.name)
        return handle.name

    def test_reads_group_ranges(self):
        """The first line is a header; each later line is name/start/stop."""
        path = self._write(
            "group start stop\n"
            "Group1 1 100\n"
            "Group2 101 250\n"
        )

        sg_dict = ts.load_group_index(path)

        self.assertEqual(sg_dict, {"Group1": [1, 100], "Group2": [101, 250]})

    def test_ignores_blank_trailing_lines(self):
        path = self._write("header\nGroup1 1 100\n\n")

        self.assertEqual(ts.load_group_index(path), {"Group1": [1, 100]})

    def test_missing_file_returns_none(self):
        """Absent index means "no sub-grouping", not a crash."""
        missing = os.path.join(
            tempfile.gettempdir(), "no_such_group_index.txt")
        self.assertFalse(os.path.exists(missing))

        self.assertIsNone(ts.load_group_index(missing))


if __name__ == '__main__':
    unittest.main()  # pragma: no cover
