import os
import unittest

from total_scattering.file_handling.load import load, create_absorption_wksp
from tests import EXAMPLE_DIR, TEST_DATA_DIR

from mantid.simpleapi import mtd


class TestLoad(unittest.TestCase):

    def setUp(self):
        self.polaris_align_and_focus_args = {
            'CalFilename': os.path.join(EXAMPLE_DIR, 'isis', 'polaris_grouping.cal'),
            'ResampleX': -6000,
            'DSpacing': False,
            'PreserveEvents': False,
            'MaxChunkSize': 8,
            'ReductionProperties': '__powderreduction'
        }

        self.nomad_align_and_focus_args = {
            'TMin': 300.0,
            'TMax': 16667.0,
            'CalFilename': os.path.join(EXAMPLE_DIR, 'sns', 'nomad_cal.h5'),
            'ResampleX': -6000,
            'DSpacing': False,
            'PreserveEvents': False,
            'MaxChunkSize': 8,
            'ReductionProperties': '__powderreduction'}

        # Highly cropped version of the workspace to improve run time
        self.si_polaris_file_path = os.path.join(TEST_DATA_DIR, 'POLARIS00097947-min.nxs')
        self.lab6_nomad_file_path = os.path.join(TEST_DATA_DIR, 'NOM_144992.nxs')

    def test_basic_load(self):
        ws_name = 'test-sample'
        actual = load(ws_name, self.si_polaris_file_path , **self.polaris_align_and_focus_args)
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
                      input_files=self.si_polaris_file_path ,
                      geometry=geometry,
                      chemical_formula=formula,
                      mass_density=mass_density,
                      **self.polaris_align_and_focus_args)
        actual = mtd[actual]
        self.assertEqual(actual.sample().getMaterial().name(), 'Si')
        mtd.clear()

    def test_load_with_abscorr(self):
        ws_name = 'test-sample'
        geometry = {'Shape': 'Cylinder',
                    'Radius': 0.3175,
                    'Center': [0.0, 0.0, 0.0],
                    'Height': 4.0}
        formula = 'Si'
        mass_density = 2.328

        a_sample, a_container = create_absorption_wksp(self.si_polaris_file_path , "SampleOnly",
                                                       geometry=geometry,
                                                       material={"ChemicalFormula": formula,
                                                                 "SampleMassDensity": mass_density})
        self.assertIsNotNone(a_sample)

        actual = load(ws_name=ws_name,
                      input_files=self.si_polaris_file_path ,
                      geometry=geometry,
                      chemical_formula=formula,
                      mass_density=mass_density,
                      absorption_wksp=a_sample,
                      **self.polaris_align_and_focus_args)
        actual = mtd[actual]
        self.assertEqual(actual.sample().getMaterial().name(), 'Si')
        mtd.clear()

    def test_load_type(self):
        ws_name = 'test-sample'
        geometry = {'Shape': 'Cylinder',
                    'Radius': 0.3,
                    'Height': 1.8}
        formula = 'La1 B6'
        mass_density = 4.72
        environment = {'Name': 'InAir', 'Container': "PAC06"}

        # Test SampleOnly
        a_sample, a_container = create_absorption_wksp(self.lab6_nomad_file_path, "SampleOnly",
                                                       geometry=geometry,
                                                       material={"ChemicalFormula": formula,
                                                                 "SampleMassDensity": mass_density},
                                                       environment=environment)
        self.assertIsNotNone(a_sample)

        actual = load(ws_name=ws_name,
                      input_files=self.lab6_nomad_file_path,
                      geometry=geometry,
                      chemical_formula=formula,
                      mass_density=mass_density,
                      absorption_wksp=a_sample,
                      **self.nomad_align_and_focus_args)
        actual = mtd[actual]
        self.assertEqual(actual.sample().getMaterial().name(), 'La1 B6')
        mtd.clear()
        
        # Test SampleAndContainer
        a_sample, a_container = create_absorption_wksp(self.lab6_nomad_file_path, "SampleAndContainer",
                                                       geometry=geometry,
                                                       material={"ChemicalFormula": formula,
                                                                 "SampleMassDensity": mass_density},
                                                       environment=environment)
        self.assertIsNotNone(a_sample)

        actual = load(ws_name=ws_name,
                      input_files=self.lab6_nomad_file_path,
                      geometry=geometry,
                      chemical_formula=formula,
                      mass_density=mass_density,
                      absorption_wksp=a_sample,
                      **self.nomad_align_and_focus_args)
        actual = mtd[actual]
        self.assertEqual(actual.sample().getMaterial().name(), 'La1 B6')
        mtd.clear()

        # Test FullPaalmanPings
        a_sample, a_container = create_absorption_wksp(self.lab6_nomad_file_path, "FullPaalmanPings",
                                                       geometry=geometry,
                                                       material={"ChemicalFormula": formula,
                                                                 "SampleMassDensity": mass_density},
                                                       environment=environment)
        self.assertIsNotNone(a_sample)

        actual = load(ws_name=ws_name,
                      input_files=self.lab6_nomad_file_path,
                      geometry=geometry,
                      chemical_formula=formula,
                      mass_density=mass_density,
                      absorption_wksp=a_sample,
                      **self.nomad_align_and_focus_args)
        actual = mtd[actual]
        self.assertEqual(actual.sample().getMaterial().name(), 'La1 B6')
        mtd.clear()

if __name__ == '__main__':
    unittest.main()  # pragma: no cover
