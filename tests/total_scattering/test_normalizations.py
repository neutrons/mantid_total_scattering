# package imports
import total_scattering.reduction.total_scattering_reduction as ts
from total_scattering.reduction.normalizations import (
    Material, calculate_and_apply_fitted_levels, to_absolute_scale, to_f_of_q)
from tests import TEST_DATA_DIR

# 3rd-party imports
from mantid.simpleapi import CloneWorkspace, DeleteWorkspace, LoadNexus
from mantid.kernel import Material as MantidMaterial
import numpy as np

# standard imports
from pathlib import Path
import unittest


class TestMaterial(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        _file_path = str(Path(TEST_DATA_DIR) / 'ceriaDP375_SofQ.nxs')
        cls.s_of_q = LoadNexus(Filename=_file_path, OutputWorkspace='s_of_q')

    @classmethod
    def tearDownClass(cls) -> None:
        DeleteWorkspace(cls.s_of_q)

    def test_init(self):
        m = Material(self.s_of_q)
        assert isinstance(m._material, MantidMaterial)

    def test_wrapper(self):
        r"""Check method calls are passed on to self._material"""
        m = Material(self.s_of_q)
        self.assertAlmostEqual(m.cohScatterLength(), 5.627, places=3)
        self.assertAlmostEqual(m.totalScatterLengthSqrd(), 32.144, places=3)
        self.assertAlmostEqual(m.totalScatterXSection(), 4.039, places=3)

    def test_properties(self):
        m = Material(self.s_of_q)
        self.assertAlmostEqual(m.bcoh_avg_sqrd, 0.31669, places=3)
        self.assertAlmostEqual(m.btot_sqrd_avg, 0.32144, places=3)
        self.assertAlmostEqual(m.laue_monotonic_diffuse_scat, 1.015, places=3)


class TestNormalizations(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        _file_path = str(Path(TEST_DATA_DIR) / 'ceriaDP375_SofQ.nxs')
        cls._s_of_q = LoadNexus(Filename=_file_path, OutputWorkspace='_s_of_q')

    @classmethod
    def tearDownClass(cls) -> None:
        DeleteWorkspace(cls._s_of_q)

    def setUp(self):
        self.s_of_q = CloneWorkspace(InputWorkspace=self._s_of_q,
                                     OutputWorkspace='s_of_q')

    def tearDown(self) -> None:
        DeleteWorkspace(self.s_of_q)

    def test_to_absolute_scale(self):
        s_of_q_norm = to_absolute_scale(self.s_of_q, 's_of_q_norm')
        m = Material(self.s_of_q)
        a, b = 1. / m.bcoh_avg_sqrd, 1 - m.laue_monotonic_diffuse_scat
        y0, y = self.s_of_q.readY(3)[7342], s_of_q_norm.readY(3)[7342]
        self.assertAlmostEqual(y, a * y0 + b, places=6)
        DeleteWorkspace('s_of_q_norm')

    def test_calculate_bank_offsets(self):
        """ Test that the bank offsets are calculated properly
        """
        config = {"SelfScatteringLevelCorrection": {"Bank2": [18.0, 20.0],
                                                    "Bank3": [20.0, 30.0],
                                                    "Bank4": [30.0, 40.0],
                                                    "Bank5": [30.0, 40.0]}}
        q_ranges = ts.get_self_scattering_level(config, 45.0)
        s_q_norm, bad_fits = calculate_and_apply_fitted_levels(self.s_of_q,
                                                               q_ranges)
        # bank 2 will have a negative offset for the given range
        self.assertEqual(len(bad_fits), 1)
        self.assertIn(2, bad_fits)
        self.assertAlmostEqual(bad_fits[2], -0.4079, places=3)
        # since bank 2 was negative, it won't get scaled (same as offset=1.0)
        offsets = {2: 1.0, 3: 0.58826, 4: 0.74313, 5: 0.79304}
        for key in q_ranges:
            norm_y = s_q_norm.readY(key - 1)
            bank_y = self.s_of_q.readY(key - 1)
            self.assertTrue(np.allclose(norm_y * offsets[key], bank_y,
                                        rtol=1e-3, equal_nan=True))

    def test_to_f_of_q(self):
        f_of_q = to_f_of_q(self.s_of_q, 'f_of_q')
        y0, y = self.s_of_q.readY(3)[7342], f_of_q.readY(3)[7342]
        c = Material(self.s_of_q).bcoh_avg_sqrd
        self.assertAlmostEqual(y, c * (y0 - 1), places=6)
        DeleteWorkspace('f_of_q')


if __name__ == '__main__':
    unittest.main()  # pragma: no cover
