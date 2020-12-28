import os
import unittest
import numpy as np

from total_scattering.file_handling.load import load, create_absorption_wksp
from tests import EXAMPLE_DIR, TEST_DATA_DIR

from mantid.simpleapi import mtd

# Expected number densities & packing fraction for tests
LAB6_NUMBER_DENSITY = 0.09764445211504061
LAB6_NUMBER_DENSITY_EFFECTIVE = 0.09764445211504061
LAB6_PACKING_FRACTION = 1.0


class TestLoad(unittest.TestCase):

    def setUp(self):
        # Set args, cropped workspace for basic load tests
        self.polaris_align_and_focus_args = {
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
        self.si_polaris_file_path = os.path.join(
            TEST_DATA_DIR,
            'POLARIS00097947-min.nxs')

        # Set args, cropped workspace, etc. for testing different abscorr types
        self.lab6_nomad_file_path = os.path.join(
            TEST_DATA_DIR,
            'NOM_144992.nxs')
        self.nomad_align_and_focus_args = {
            'CalFilename': os.path.join(EXAMPLE_DIR, 'sns', 'nomad_cal.h5'),
            'ResampleX': -6000,
            'DSpacing': False,
            'PreserveEvents': False,
            'MaxChunkSize': 8,
            'TMin': 300.0,
            'TMax': 16667.0,
            'ReductionProperties': '__powderreduction'}

        self.type_test_geometry = {
            'Shape': 'Cylinder',
            'Radius': 0.3,  # measured in cm
            'Height': 1.8}  # measured in cm
        self.type_test_environment = {
            'Name': 'InAir',
            'Container': "PAC06"}
        self.type_test_material = {
            'ChemicalFormula': 'La1 B6',
            'SampleMassDensity': 4.72}

    def test_basic_load(self):
        ws_name = 'test-sample'
        actual = load(
            ws_name,
            self.si_polaris_file_path,
            **self.polaris_align_and_focus_args)
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
        actual = load(
            ws_name=ws_name,
            input_files=self.si_polaris_file_path,
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

        a_sample, a_container = create_absorption_wksp(
            self.si_polaris_file_path,
            "SampleOnly",
            geometry=geometry,
            material={
                "ChemicalFormula": formula,
                "SampleMassDensity": mass_density})

        self.assertIsNotNone(a_sample)

        actual = load(
            ws_name=ws_name,
            input_files=self.si_polaris_file_path,
            geometry=geometry,
            chemical_formula=formula,
            mass_density=mass_density,
            absorption_wksp=a_sample,
            **self.polaris_align_and_focus_args)

        actual = mtd[actual]
        self.assertEqual(actual.sample().getMaterial().name(), 'Si')
        mtd.clear()

    def test_load_sampleonly(self):
        ws_name = 'test-sample'

        a_sample, a_container = create_absorption_wksp(
            self.lab6_nomad_file_path,
            "SampleOnly",
            geometry=self.type_test_geometry,
            material=self.type_test_material,
            environment=self.type_test_environment)

        self.assertIsNotNone(a_sample)

        actual = load(
            ws_name=ws_name,
            input_files=self.lab6_nomad_file_path,
            geometry=self.type_test_geometry,
            chemical_formula=self.type_test_material['ChemicalFormula'],
            mass_density=self.type_test_material['SampleMassDensity'],
            absorption_wksp=a_sample,
            **self.nomad_align_and_focus_args)

        actual = mtd[actual]

        volume = actual.sample().getShape().volume()
        assert volume == np.pi * np.square(0.003) * 0.018

        material = actual.sample().getMaterial()
        assert material.name() == 'La1 B6'
        assert material.numberDensity == LAB6_NUMBER_DENSITY
        assert material.numberDensityEffective == LAB6_NUMBER_DENSITY_EFFECTIVE
        assert material.packingFraction == LAB6_PACKING_FRACTION
        mtd.clear()

    def test_load_samplecontainer(self):
        ws_name = 'test-sample'

        a_sample, a_container = create_absorption_wksp(
            self.lab6_nomad_file_path,
            "SampleAndContainer",
            geometry=self.type_test_geometry,
            material=self.type_test_material,
            environment=self.type_test_environment)

        self.assertIsNotNone(a_sample)

        actual = load(
            ws_name=ws_name,
            input_files=self.lab6_nomad_file_path,
            geometry=self.type_test_geometry,
            chemical_formula=self.type_test_material['ChemicalFormula'],
            mass_density=self.type_test_material['SampleMassDensity'],
            absorption_wksp=a_sample,
            **self.nomad_align_and_focus_args)

        actual = mtd[actual]

        volume = actual.sample().getShape().volume()
        assert volume == np.pi * np.square(0.003) * 0.018

        material = actual.sample().getMaterial()
        assert material.name() == 'La1 B6'
        assert material.numberDensity == LAB6_NUMBER_DENSITY
        assert material.numberDensityEffective == LAB6_NUMBER_DENSITY_EFFECTIVE
        assert material.packingFraction == LAB6_PACKING_FRACTION
        mtd.clear()

    def test_load_fullpaalmanpings(self):
        ws_name = 'test-sample'

        a_sample, a_container = create_absorption_wksp(
            self.lab6_nomad_file_path,
            "FullPaalmanPings",
            geometry=self.type_test_geometry,
            material=self.type_test_material,
            environment=self.type_test_environment)

        self.assertIsNotNone(a_sample)

        actual = load(
            ws_name=ws_name,
            input_files=self.lab6_nomad_file_path,
            geometry=self.type_test_geometry,
            chemical_formula=self.type_test_material['ChemicalFormula'],
            mass_density=self.type_test_material['SampleMassDensity'],
            absorption_wksp=a_sample,
            **self.nomad_align_and_focus_args)

        actual = mtd[actual]

        volume = actual.sample().getShape().volume()
        assert volume == np.pi * np.square(0.003) * 0.018

        material = actual.sample().getMaterial()
        assert material.name() == 'La1 B6'
        assert material.numberDensity == LAB6_NUMBER_DENSITY
        assert material.numberDensityEffective == LAB6_NUMBER_DENSITY_EFFECTIVE
        assert material.packingFraction == LAB6_PACKING_FRACTION
        mtd.clear()


if __name__ == '__main__':
    unittest.main()  # pragma: no cover
