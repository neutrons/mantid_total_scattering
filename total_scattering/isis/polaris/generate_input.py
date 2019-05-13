import os


from total_scattering.utils import ROOT_DIR
POLARIS_DIR = os.path.join(ROOT_DIR, 'total_scattering', 'isis', 'polaris')


JSON_MSG = """
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
}
"""


def generate_input_json():
    calib_path = os.path.join(ROOT_DIR, 'examples', 'isis', 'polaris_grouping.cal')
    character_path = os.path.join(POLARIS_DIR, 'character.txt')
    cache_path = os.path.join(ROOT_DIR, 'cache')
    output_path = os.path.join(ROOT_DIR, 'output')
    if os.name == 'nt':
        calib_path = calib_path.replace('\\', '\\\\')
        character_path = character_path.replace('\\', '\\\\')
        cache_path = cache_path.replace('\\', '\\\\')
        output_path = output_path.replace('\\', '\\\\')
    with open(os.path.join(POLARIS_DIR, 'test_input.json'), 'w') as input_file:
        input_file.write(JSON_MSG % (calib_path,
                                     character_path,
                                     cache_path,  # cache dir
                                     output_path))  # output dir


def clean_up():
    os.removedirs(os.path.join(ROOT_DIR, 'cache'))
    os.removedirs(os.path.join(ROOT_DIR, 'output'))
    os.remove(os.path.join(POLARIS_DIR, 'test_input.json'))
