import os
import unittest
import numpy as np

from total_scattering.file_handling.load import load
from total_scattering.file_handling.save import save_banks, save_file
from tests import EXAMPLE_DIR, TEST_DATA_DIR

from mantid.simpleapi import mtd, \
    LoadNexusProcessed, LoadAscii, ConvertToHistogram


class TestSave(unittest.TestCase):

    def setUp(self):
        align_and_focus_args = {
            'CalFilename': os.path.join(EXAMPLE_DIR, 'isis', 'polaris_grouping.cal'),
            'ResampleX': -6000,
            'DSpacing': False,
            'PreserveEvents': False,
            'MaxChunkSize': 8,
            'ReductionProperties': '__powderreduction'
        }
        # Highly cropped version of the workspace to improve run time
        ws_name = 'test-sample'
        sample_file_path = os.path.join(TEST_DATA_DIR, 'POLARIS00097947-min.nxs')
        wksp = load(ws_name, sample_file_path, **align_and_focus_args)

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
        self.assertTrue(os.path.isfile(self.out_nxs))
        mtd.clear()

    def test_save_banks_relative_path(self):
        save_banks(self.wksp, self.out_nxs, 'wksp', './output')
        self.assertTrue(os.path.isfile(os.path.join('./output', self.out_nxs)))
        mtd.clear()

    def test_save_banks_check_contents(self):
        save_banks(self.wksp, self.out_nxs, 'wksp', '.')
        out_wksp = LoadNexusProcessed(self.out_nxs)
        self.assertEqual(out_wksp.blocksize(),
                         self.wksp.blocksize())
        self.assertEqual(out_wksp.getNumberHistograms(),
                         self.wksp.getNumberHistograms())
        self.assertTrue(np.array_equal(out_wksp.getAxis(0).extractValues(),
                                       self.wksp.getAxis(0).extractValues())
                        )

    def test_save_banks_binning(self):
        save_banks(self.wksp, self.out_nxs, 'wksp', '.', Binning='0,100,10000')
        out_wksp = LoadNexusProcessed(self.out_nxs)
        self.assertNotEqual(out_wksp.blocksize(),
                            self.wksp.blocksize())
        self.assertEqual(out_wksp.blocksize(), 100)

    def test_save_banks_grouping(self):
        # TODO: Will have to implement when we have event test data.
        '''
        Below does not work since this POLARIS test data is Histogrammed
        and does not contain counts

        will need to re-add CreateGroupingWorkspace to imported algorithms

        # Create grouping workspace.
        grp_ws, nspectra, grp_count = CreateGroupingWorkspace(InstrumentName='POLARIS',
                                                              GroupDetectorsBy="All")
        save_banks(self.wksp, self.out_nxs, 'wksp_title', '.', GroupingWorkspace=grp_ws)
        out_wksp = LoadNexusProcessed(self.out_nxs)

        self.assertEqual(grp_ws.blocksize(), 1)
        self.assertEqual(nspectra, 3008)
        self.assertEqual(grp_count, 1)
        self.assertNotEqual(out_wksp.getNumberHistograms(),
                            self.wksp.getNumberHistograms())
        self.assertEqual(out_wksp.getNumberHistograms(), 1)
        '''

    def test_save_file_exists(self):
        save_file(self.wksp, self.out_ascii)
        self.assertTrue(os.path.isfile(self.out_ascii))

    def test_save_file_check_contents(self):
        save_file(self.wksp, self.out_ascii)
        out_wksp = LoadAscii(self.out_ascii, Separator='Space')
        out_wksp = ConvertToHistogram(out_wksp)

        self.assertEqual(out_wksp.blocksize(),
                         self.wksp.blocksize())
        self.assertEqual(out_wksp.getNumberHistograms(),
                         self.wksp.getNumberHistograms())
        self.assertTrue(np.allclose(out_wksp.getAxis(0).extractValues(),
                                    self.wksp.getAxis(0).extractValues())
                        )


if __name__ == '__main__':
    unittest.main()  # pragma: no cover
