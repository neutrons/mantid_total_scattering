import json
import os

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

# Configurations shipped with the package, used when the facility location
# above is not reachable (container, laptop, Galaxy job, ...).
PACKAGED_CONFIG_DIR = os.path.join(
    os.path.abspath(os.path.dirname(__file__)), "config")

# Environment overrides, checked before the facility location.
CONFIG_FILE_ENV_VAR = "MTS_AUTO_CONFIG_FILE"
CONFIG_DIR_ENV_VAR = "MTS_CONFIG_DIR"

CONFIG_BASENAME = "auto_config.json"


def _packaged_config_candidates(facility, instrument):
    """Packaged configs, most specific first: instrument then facility."""
    return [
        os.path.join(PACKAGED_CONFIG_DIR, facility, instrument,
                     CONFIG_BASENAME),
        os.path.join(PACKAGED_CONFIG_DIR, facility, CONFIG_BASENAME),
    ]


def resolve_config_file(facility, instrument, config_loc_in=None,
                        config_file=None):
    """Find the facility/instrument configuration file.

    The candidates are tried in order and the first one that exists wins:

    1. ``config_file``, passed in explicitly (e.g. from the input JSON)
    2. the file named by ``$MTS_AUTO_CONFIG_FILE``
    3. ``auto_config.json`` under ``$MTS_CONFIG_DIR``
    4. the facility location (the SNS analysis cluster share)
    5. the copy packaged with ``total_scattering``

    @param facility: str, e.g. "SNS"
    @param instrument: str, e.g. "NOM"
    @param config_loc_in: dict, overrides the module level ``config_loc``
    @param config_file: str, an explicit path to the configuration file

    @return: str, path to an existing configuration file
    """
    if config_loc_in is None:
        config_loc_in = config_loc

    candidates = []
    if config_file:
        candidates.append(config_file)

    env_file = os.environ.get(CONFIG_FILE_ENV_VAR)
    if env_file:
        candidates.append(env_file)

    env_dir = os.environ.get(CONFIG_DIR_ENV_VAR)
    if env_dir:
        candidates.append(os.path.join(env_dir, CONFIG_BASENAME))

    facility_file = config_loc_in.get(facility, {}).get(instrument)
    if facility_file:
        candidates.append(facility_file)

    candidates.extend(_packaged_config_candidates(facility, instrument))

    for candidate in candidates:
        if os.path.isfile(candidate):
            return candidate

    raise FileNotFoundError(
        "Could not find a configuration file for facility '{0}' and "
        "instrument '{1}'. Looked in: {2}. Set '{3}' or '{4}' to point at "
        "one, or pass 'AutoConfigFile' in the reduction input.".format(
            facility, instrument, ", ".join(candidates),
            CONFIG_FILE_ENV_VAR, CONFIG_DIR_ENV_VAR))


class ParamsLoader:
    def __init__(self, facility, instrument, config_loc_in=None,
                 config_file=None):
        if config_loc_in is None:
            config_loc_in = config_loc

        self.config_loc = config_loc_in
        self.abs_ms_sn = abs_ms_sn

        self.config_file = resolve_config_file(
            facility, instrument,
            config_loc_in=config_loc_in,
            config_file=config_file)
        # Several auxiliary files (grouping, group index) are looked up
        # alongside the configuration file, so the directory it came from
        # matters as much as its contents.
        self.config_dir = os.path.dirname(self.config_file)

        with open(self.config_file) as f:
            self.config_params = json.load(f)
