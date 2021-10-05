# standard imports
import unittest
from total_scattering.reduction.validator import validateConfig


class CliTest(unittest.TestCase):
    demoConfig = {
        "Facility": "SNS",
        "Instrument": "NOM",
        "Title": "ceriaDP375_PAC06_at_300K_IDL_Calib",
        "Calibration": {
            "Filename": "test.h5"  # don't need a real file here
        },
        "AlignAndFocusArgs": {"TMin": 300.0, "TMax": 16667.0},
        "HighQLinearFitRange": 0.60,
        "Merging": {"QBinning": [0.0, 0.005, 40.0], "SumBanks": [3]},
        "Sample": {},
        "Normalization": {},
        "CacheDir": "./tmp",
        "OutputDir": "./output",
    }

    def test_validateConfig(self):
        # case 0: no abs nor ms correction
        validateConfig(self.demoConfig)

        # case 1: abs correction, but no ms correction for sample
        config = self.demoConfig.copy()
        config["Sample"] = {"AbsorptionCorrection": {"Type": "SampleOnly"}}
        validateConfig(config)

        # case 2: abs and ms correction for sample only
        config = self.demoConfig.copy()
        config["Sample"] = {
            "AbsorptionCorrection": {"Type": "SampleOnly"},
            "MultipleScatteringCorrection": {"Type": "SampleOnly"},
        }
        validateConfig(config)

        # case 3: abs and ms correction for sample and container
        config = self.demoConfig.copy()
        config["Sample"] = {
            "AbsorptionCorrection": {"Type": "FullPaalmanPings"},
            "MultipleScatteringCorrection": {"Type": "SampleAndContainer"},
        }
        config["Normalization"] = {
            "AbsorptionCorrection": {"Type": "SampleOnly"},
            "MultipleScatteringCorrection": {"Type": "SampleOnly"},
        }


if __name__ == "__main__":
    unittest.main()
