import os
from os.path import abspath, dirname, join

ROOT_DIR = abspath(join(dirname(abspath(__file__)), '..'))
EXAMPLE_DIR = os.path.join(ROOT_DIR, 'examples')
TEST_DATA_DIR = os.path.join(ROOT_DIR, 'tests', 'data')
