import unittest
import total_scattering.reduction.total_scattering_reduction as ts


class TestUtilsForReduction(unittest.TestCase):

    def setUp(self):
        pass

    def test_one_and_only_one(self):
        inputList = [True, False, False]
        output = ts.one_and_only_one(inputList)
        self.assertTrue(output)