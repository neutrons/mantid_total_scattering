# standard imports
import random
import string
from os.path import abspath, dirname, join

ROOT_DIR = abspath(join(dirname(abspath(__file__)), '..'))


def random_name(n: int = 12, prefix: str = '_', suffix: str = '') -> str:
    r"""
    Random word intended as the name of a temporary workspace

    Parameters
    ----------
    n: length of the word, not counting the prefix and suffix
    prefix: prepend this string to the random word
    suffix: append this string to the random word
    """
    letters = [random.choice(string.ascii_letters) for _ in range(n)]  # nosec
    return prefix + ''.join(letters) + suffix
