r"""
Functions assisting in tests writing
"""
# standard imports
import random
import string
from typing import Any, Dict
import unittest


def compare(expected: Dict[str, Any],
            actual: Dict[str, Any]) -> None:
    r"""Assert values of two dictionaries

    :param expected: dictionary containing the correct values
    :param actual: dictionary containing values compared against correct ones
    """
    assertor = unittest.TestCase()
    comparators = {
        str: assertor.assertEqual,
        int: assertor.assertEqual,
        float: assertor.assertAlmostEqual
    }
    for key, value in expected.items():
        if isinstance(value, dict) and 'Kwargs' in value:
            comparator = comparators[type(value['Value'])]
            comparator(value['Value'], actual[key], **value['Kwargs'])
        else:
            comparator = comparators[type(value)]
            comparator(value, actual[key])
