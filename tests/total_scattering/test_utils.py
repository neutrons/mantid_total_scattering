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
            ts.extract_key_match_from_dict(keys, dictionary)

    def test_get_sample_when_match(self):
        """ Test extracting sample info from config
        """
        config = {"Sample": {"Runs": "10-20"}}
        sample_value = ts.get_sample(config)
        self.assertEqual(config["Sample"], sample_value)

    def test_get_sample_raise_error_when_no_match(self):
        """ Test that we raise an error when no Sample key found
        """
        config = {"BadKey": {"Runs": "10-20"}}
        with self.assertRaises(Exception):
            ts.get_sample(config)

    def test_get_normalization_when_match_for_vanadium(self):
        """ Test extracting vanadium info from config
        """
        config = {"Vanadium": {"Runs": "10-20"}}
        norm_value = ts.get_normalization(config)
        self.assertEqual(config["Vanadium"], norm_value)

    def test_get_normalization_when_match_for_normalization(self):
        """ Test extracting normalization info from config
        """
        config = {"Normalization": {"Runs": "10-20"}}
        norm_value = ts.get_normalization(config)
        self.assertEqual(config["Normalization"], norm_value)

    def test_get_normalization_when_match_for_normalisation(self):
        """ Test extracting normalisation info from config
        """
        config = {"Normalisation": {"Runs": "10-20"}}
        norm_value = ts.get_normalization(config)
        self.assertEqual(config["Normalisation"], norm_value)

    def test_get_normalization_raise_error_when_no_match(self):
        """ Test that we raise an error when no normalization keys found
        """
        config = {"BadKey": {"Runs": "10-20"}}
        with self.assertRaises(Exception):
            ts.get_normalization(config)
