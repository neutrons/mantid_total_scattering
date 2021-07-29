import total_scattering.reduction.total_scattering_reduction as ts
from total_scattering.reduction.normalizations import calculate_fitted_levels
from tests import TEST_DATA_DIR
from mantid.simpleapi import Load

import os
import numpy as np
import unittest


class TestNormalizations(unittest.TestCase):

    def setUp(self):
        pass

    def test_calculate_bank_offsets(self):
        """ Test that the bank offsets are calculated properly
        """
        s_q_wksp = Load(Filename=os.path.join(TEST_DATA_DIR,
                                              "ceriaDP375_SofQ.nxs"))
        config = {"SelfScatteringLevelCorrection": {"Bank3": [20.0, 30.0],
                                                    "Bank4": [30.0, 40.0],
                                                    "Bank5": [30.0, 40.0]}}
        q_ranges = ts.get_self_scattering_level(config, 45.0)
        s_q_norm = calculate_fitted_levels(s_q_wksp, q_ranges)
        offsets = {3: 0.58826, 4: 0.74313, 5: 0.79304}
        for key in q_ranges:
            norm_y = s_q_norm.readY(key - 1)
            bank_y = s_q_wksp.readY(key - 1)
            self.assertTrue(np.allclose(norm_y * offsets[key], bank_y,
                                        rtol=1e-3, equal_nan=True))


if __name__ == '__main__':
    unittest.main()  # pragma: no cover
