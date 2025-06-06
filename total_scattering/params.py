import json

config_loc = {
    "SNS": {
        "NOM": "/SNS/NOM/shared/autoreduce/configs/auto_config.json",
        "PG3": "/SNS/PG3/shared/autoreduce/configs/auto_config.json",
        "CORELLI": "/SNS/NOM/shared/config/CORELLI/auto_config_corelli.json"
    }
}

abs_ms_sn = {
    "SampleOnly": "SO",
    "SampleAndContainer": "SC",
    "FullPaalmanPings": "FPP"
}


class ParamsLoader:
    def __init__(self, facility, instrument, config_loc_in=None):
        if config_loc_in is None:
            config_loc_in = config_loc

        self.config_loc = config_loc_in
        self.abs_ms_sn = abs_ms_sn

        with open(self.config_loc[facility][instrument]) as f:
            self.config_params = json.load(f)
