import json
import os
import unittest

import mantid_total_scattering as total_scattering
from utils import ROOT_DIR


class TestPolarisTotalScattering(unittest.TestCase):

    def test_silicon(self):
        """
        Run polaris silicon data through total scattering script
        """
        with open(os.path.join(ROOT_DIR, 'isis', 'polaris', 'input.json'), 'r') as handle:
            config = json.load(handle)
        actual = total_scattering.main(config)
        self.assertEqual(actual.getNumberHistograms(), 5)
        self.assertEqual(actual.getTitle(), "10: Silicon 640b, 700MeV, chopper stopped")
        self.assertEqual(actual.getInstrument().getFullName(), "POLARIS")