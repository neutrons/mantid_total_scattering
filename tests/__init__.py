import os
from total_scattering.utils import ROOT_DIR

# https://docs.travis-ci.com/user/environment-variables/#default-environment-variables
IN_TRAVIS = os.getenv('TRAVIS', False)
EXAMPLE_DIR = os.path.join(ROOT_DIR, 'examples')
TEST_DATA_DIR = os.path.join(ROOT_DIR, 'tests', 'data')
