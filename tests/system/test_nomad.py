r"""

"""

# package imports
from total_scattering.reduction import TotalScatteringReduction
from tests import EXAMPLE_DIR
from tests.utils import compare

# standard imports
import glob
import json
from pathlib import Path
import unittest


class NomadTotalScatteringSystemTest(unittest.TestCase):

    NOM_DIR_ACCESIBLE = Path('/SNS/NOM').exists()

    @unittest.skipIf(not NOM_DIR_ACCESIBLE, 'Do not run system test on build servers')
    def test_examples(self):
        pattern = str(Path(EXAMPLE_DIR) / 'sns' / 'nomad_*.json')
        for file_json in glob.iglob(pattern):

            # reduce only when a file with expected values exists
            file_expected = file_json.replace('.json', '.expt')  # JSON file containing expected values
            if not Path(file_expected).exists():
                continue

            # reduce the file
            with open(file_json, 'r') as handle:
                config = json.load(handle)
                workspace = TotalScatteringReduction(config)  # Workspace

            # create dictionary of actual values from the workspace
            actual = {
                'title': workspace.getTitle(),
                'instrument': workspace.getInstrument().getFullName()
                # TODO add more meaninful comparisons
            }

            # compare actual with expected values
            with open(file_expected, 'r') as handle:
                expected = json.load(handle)  # dictionary of expected values
            compare(expected, actual)


if __name__ == '__main__':
    unittest.main()  # pragma: no cover