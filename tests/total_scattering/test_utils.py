import unittest
import total_scattering.reduction.total_scattering_reduction as ts


class TestUtilsForReduction(unittest.TestCase):

    def setUp(self):
        pass

    def test_one_and_only_one(self):
        """ Test for the one and only one true value utility function
        """
        inputList = [True, False, False]
        output = ts.one_and_only_one(inputList)
        self.assertTrue(output)

    def test_find_key_match_in_dict(self):
        """ Test for using list of keys to get back value of
        one and only one matching key from dict
        """
        keys = ["Vanadium", "Normalization", "Normalisation"]
        dictionary = {"Normalization": True}
        match_value = ts.find_key_match_in_dict(keys, dictionary)
        self.assertTrue(match_value)

    def test_extract_key_match_from_dict_for_match(self):
        """ Test for using list of keys to get back value of
        one and only one matching key from dict with error handling wrapper
        """
        keys = ["Vanadium", "Normalization", "Normalisation"]
        dictionary = {"Normalization": True}
        match_value = ts.extract_key_match_from_dict(keys, dictionary)
        self.assertTrue(match_value)

    def test_extract_key_match_from_dict_raise_error_when_no_match(self):
        """ Test that we raise an error when no keys found in dict
        """
        keys = ["Vanadium", "Normalization", "Normalisation"]
        dictionary = {"Sample": True}
        with self.assertRaises(Exception):
            match_value = ts.extract_key_match_from_dict(keys, dictionary)
        