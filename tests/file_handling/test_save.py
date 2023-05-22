import os
import unittest

from total_scattering.file_handling.load import load
from total_scattering.file_handling.save import save_banks, save_file
from tests import EXAMPLE_DIR, TEST_DATA_DIR
from mantid.simpleapi import mtd


class TestSave(unittest.TestCase):

    def setUp(self):
        align_and_focus_args = {
            'CalFilename': os.path.join(
                EXAMPLE_DIR,
                'isis',
                'polaris_grouping.cal'),
            'ResampleX': -6000,
            'DSpacing': False,
            'PreserveEvents': False,
            'MaxChunkSize': 8,
            'ReductionProperties': '__powderreduction'
        }
        # Highly cropped version of the workspace to improve run time
        ws_name = 'test-sample'
        sample_file_path = os.path.join(
            TEST_DATA_DIR,
            'POLARIS00097947-min.nxs')
        wksp = load(ws_name, sample_file_path, None, **align_and_focus_args)

        self.wksp = mtd[wksp]
        self.out_nxs = '%s.nxs' % ws_name
        self.out_ascii = '%s.nxs' % ws_name

    def tearDown(self):
        mtd.clear()
        if os.path.isfile(self.out_nxs):
            os.remove(self.out_nxs)
        if os.path.isfile(self.out_ascii):
            os.remove(self.out_ascii)

    def test_save_banks_exists(self):
        save_banks(self.wksp, self.out_nxs, 'wksp', '.')
        self.assertTrue(os.path.isfile(os.path.join(".", "SofQ", self.out_nxs)))
        mtd.clear()

    def test_save_banks_relative_path(self):
        save_banks(self.wksp, self.out_nxs, 'wksp', './output')
        self.assertTrue(os.path.isfile(os.path.join('./output', "SofQ", self.out_nxs)))
        mtd.clear()

    def test_save_file_exists(self):
        save_file(self.wksp, self.out_ascii)
        self.assertTrue(os.path.isfile(self.out_ascii))


if __name__ == '__main__':
    unittest.main()  # pragma: no cover
