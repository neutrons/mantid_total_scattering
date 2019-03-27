import json
import os
import unittest

from total_scattering.isis.polaris.generate_input import generate_input_json
from total_scattering.utils import ROOT_DIR
from tests import IN_TRAVIS
try:
    import mantid_total_scattering as total_scattering  # TODO will not work that is why it is in here
except ImportError:
    total_scattering = None  # mark as missing

class PolarisTotalScatteringSystemTest(unittest.TestCase):
    @unittest.skipIf(IN_TRAVIS or not total_scattering, 'Do not run thest on build servers')
    def test_silicon(self):
        """
        Run polaris silicon data through total scattering script
        """
        generate_input_json()
        with open(os.path.join(ROOT_DIR, 'total_scattering', 'isis', 'polaris', 'test_input.json'), 'r') as handle:
            config = json.load(handle)
        actual = total_scattering.main(config)
        self.assertEqual(actual.getNumberHistograms(), 5)
        self.assertEqual(actual.getTitle(), "10: Silicon 640b, 700MeV, chopper stopped")
        self.assertEqual(actual.getInstrument().getFullName(), "POLARIS")

if __name__ == '__main__':
    unittest.main()  # pragma: no cover
