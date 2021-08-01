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

    def test_get_self_scattering_level(self):
        """ Test that a dictionary of self scattering corrections is properly
        returned for each bank number
        """
        config = {"SelfScatteringLevelCorrection": {"Bank3": [20.0, 30.0],
                                                    "Bank4": [30.0, 40.0],
                                                    "Bank5": [30.0, 45.0]}}
        self_scattering = ts.get_self_scattering_level(config, 45.0)
        self.assertIsInstance(self_scattering, dict)
        self.assertEqual(len(self_scattering), 3)
        self.assertIn(3, self_scattering)
        self.assertIn(4, self_scattering)
        self.assertIn(5, self_scattering)
        self.assertEqual(self_scattering[3], (20.0, 30.0))

    def test_get_self_scattering_level_invalid_bank_name(self):
        """ Test the self scattering option with an invalid bank name
        """
        config = {"SelfScatteringLevelCorrection": {"bank1": [20.0, 30.0]}}
        with self.assertRaises(RuntimeError):
            ts.get_self_scattering_level(config, 45.0)

    def test_get_self_scattering_level_min_max(self):
        """ Test the self scattering option where bank level min > max
        """
        config = {"SelfScatteringLevelCorrection": {"Bank1": [20.0, 10.0]}}
        with self.assertRaises(RuntimeError):
            ts.get_self_scattering_level(config, 45.0)

    def test_get_self_scattering_level_qbin_clamp(self):
        """ Test that the self scattering max values are clamped by the
        max q-binning value
        """
        config = {"SelfScatteringLevelCorrection": {"Bank5": [30.0, 45.0]}}
        self_scattering = ts.get_self_scattering_level(config, 40.0)
        self.assertEqual(self_scattering[5], (30.0, 40.0))


if __name__ == '__main__':
    unittest.main()  # pragma: no cover
