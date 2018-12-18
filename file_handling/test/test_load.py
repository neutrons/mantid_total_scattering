import os
import unittest

from file_handling.load import load
from utils import ROOT_DIR

from mantid.simpleapi import mtd


class TestLoad(unittest.TestCase):

    def setUp(self):
        align_and_focus_args = {
            'CalFilename': os.path.join(ROOT_DIR, 'isis', 'polaris', 'grouping.cal'),
            'ResampleX': -6000,
            'DSpacing': False,
            'PreserveEvents': False,
            'MaxChunkSize': 8,
            'ReductionProperties': '__powderreduction'
        }

        self.loadKwargs = {'AlignAndFocusArgs': align_and_focus_args}

        # Highly cropped version of the workspace to improve run time
        self.sample_file_path = os.path.join(ROOT_DIR, 'test_data', 'POLARIS00097947-min.nxs')

    def test_basic_load(self):
        ws_name = 'test-sample'
        actual = load(ws_name, self.sample_file_path, **self.loadKwargs)
        actual = mtd[actual]
        self.assertEqual(actual.name(), ws_name)
        self.assertEqual(actual.getNumberHistograms(), 5)
        self.assertEqual(actual.getAxis(0).getUnit().caption(), 'q')
        mtd.clear()

    def test_load_with_material(self):
        ws_name = 'test-sample'
        geometry = {'Shape': 'Cylinder',
                    'Radius': 0.3175,
                    'Center': [0.0, 0.0, 0.0],
                    'Height': 4.0}
        formula = 'Si'
        mass_density = 2.328
        actual = load(ws_name=ws_name,
                      input_files=self.sample_file_path,
                      geometry=geometry,
                      chemical_formula=formula,
                      mass_density=mass_density,
                      **self.loadKwargs)
        actual = mtd[actual]
        self.assertEqual(actual.sample().getMaterial().name(), 'Si')
        mtd.clear()

    # TODO: Cover the correlation chopper with test.
    # Need Corelli data uploaded to remote test server first.


if __name__ == '__main__':
    unittest.main()  # pragma: no cover
