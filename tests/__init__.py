import os
from os.path import abspath, dirname, join

# https://docs.travis-ci.com/user/environment-variables/#default-environment-variables
IN_TRAVIS = os.getenv('TRAVIS', False)
ROOT_DIR = abspath(join(dirname(abspath(__file__)), '..'))
EXAMPLE_DIR = os.path.join(ROOT_DIR, 'examples')
TEST_DATA_DIR = os.path.join(ROOT_DIR, 'tests', 'data')
