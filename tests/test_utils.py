
r"""

"""
# package imports
from tests.utils import compare

# standard imports
import unittest

class UtilsTest(unittest.TestCase):

    def test_compare(self):
        expected = {
            'one': 'one',
            'two': 2,
            'three_almost': 2.99,
            'four_almost': {'Value': 3.99, 'Kwargs': dict(places=3)}
        }
        actual = {}
        for key, value in expected.items():
            if isinstance(value, dict) and 'Kwargs' in value:
                actual[key] = value['Value']  # strip Kwargs from the value associated to the key
            else:
                actual[key] = value
        compare(expected, actual)


if __name__ == '__main__':
    unittest.main()  # pragma: no cover
