import os


from total_scattering.utils import ROOT_DIR


def generate_input_json():
    with open(os.path.join(ROOT_DIR, 'isis', 'polaris', 'test_input.json'), 'w') as input_file:
        input_file.write(
            """
            {
            "Facility": "ISIS",
            "Instrument": "POLARIS",
            "Title": "Silicon (NIST SRM 640b)",
            "Sample": {"Runs": "97947",
                       "Background": {"Runs": "97948-97951",
                                      "Background": {"Runs": "97952-97954"}
                                      },
                       "Material": "Si",
                       "MassDensity": 2.328,
                       "PackingFraction": 0.35,
                       "Geometry": {"Radius": 0.3175,
                                    "Height": 4.0},
                       "AbsorptionCorrection": {"Type": "Carpenter"},
                       "MultipleScatteringCorrection": {"Type": "Carpenter"},
                       "InelasticCorrection": {"Type": "Placzek",
                                               "Order": "1st",
                                               "Self": true,
                                               "Interference": false,
                                               "FitSpectrumWith": "GaussConvCubicSpline",
                                               "LambdaBinningForFit": "0.16,0.04,2.8",
                                               "LambdaBinningForCalc": "0.16,0.0001,2.9"}
                       },
            "Vanadium": {"Runs": "97942-97946",
                         "Background": {"Runs": "97952-97954"},
                         "Material": "V",
                         "MassDensity": 6.11,
                         "PackingFraction": 1.0,
                         "Geometry": {"Radius": 0.025,
                                      "Height": 4.0},
                         "AbsorptionCorrection": {"Type": "Carpenter"},
                         "MultipleScatteringCorrection": {"Type": "Carpenter"},
                         "InelasticCorrection": {"Type": "Placzek",
                                                 "Order": "1st",
                                                 "Self": true,
                                                 "Interference": false,
                                                 "FitSpectrumWith": "GaussConvCubicSpline",
                                                 "LambdaBinningForFit": "0.16,0.04,2.8",
                                                 "LambdaBinningForCalc": "0.1,0.0001,3.0"}
                         },
            "Calibration": {"Filename": "%s"},
            "HighQLinearFitRange": 0.60,
            "Merging": {"QBinning": [0.35, 0.05, 11.4],
                        "SumBanks": [3],
                        "Characterizations": {"Filename": "%s"}
                        },
            "CacheDir": "%s",
            "OutputDir": "%s"
            }""" % (os.path.join(ROOT_DIR, 'isis', 'polaris', 'grouping.cal'),  # Calibration grouping file
                    os.path.join(ROOT_DIR, 'isis', 'polaris', 'character.txt'),  # Characterisation file
                    ROOT_DIR,  # cache dir
                    ROOT_DIR)  # output dir
        )
