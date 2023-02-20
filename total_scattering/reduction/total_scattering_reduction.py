#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function)

import os
import itertools
import numpy as np
from scipy.constants import Avogadro

from mantid import mtd
from mantid.kernel import Logger
from mantid.simpleapi import \
    CarpenterSampleCorrection, \
    CloneWorkspace, \
    CompressEvents, \
    ConvertToDistribution, \
    ConvertToHistogram,\
    ConvertUnits, \
    CreateEmptyTableWorkspace, \
    CreateGroupingWorkspace, \
    CropWorkspaceRagged, \
    Divide, \
    FFTSmooth, \
    GenerateEventsFilter, \
    GroupWorkspaces, \
    Load, \
    LoadDetectorsGroupingFile, \
    LoadDiffCal, \
    MayersSampleCorrection, \
    Minus, \
    PDDetermineCharacterizations, \
    PropertyManagerDataService, \
    Rebin, \
    RebinToWorkspace, \
    SaveGSS, \
    SaveFocusedXYE, \
    SetSample, \
    SetUncertainties, \
    StripVanadiumPeaks, \
    CalculateEfficiencyCorrection, \
    CalculatePlaczek, \
    LoadNexusMonitors, \
    ConvertToPointData, \
    Multiply, \
    GetIPTS, \
    SaveNexus, \
    LoadNexus, \
    RenameWorkspace

from total_scattering.file_handling.load import load, create_absorption_wksp
from total_scattering.file_handling.save import save_banks
from total_scattering.inelastic.placzek import FitIncidentSpectrum
from total_scattering.reduction.normalizations import (
    Material, calculate_and_apply_fitted_levels, to_absolute_scale, to_f_of_q)


# Utilities
def generate_cropping_table(qmin, qmax):
    ''' Generate a Table workspace that can be used to crop
    another workspace in reciprocal space (ie MomentumTransfer)

    :param qmin: list of Qmin values to crop each spectra by
    :type qmin: str (example: '0.2,0.4,1.0')
    :param qmax: list of Qmax values to crop each spectra by
    :type qmax: str (example: '10.5,12.0,40.0')

    :return: Cropping table with columns of ("SpectraList","Xmin","Xmax")
    :rtype: TableWorkspace
    '''
    mask_info = CreateEmptyTableWorkspace()
    mask_info.addColumn("str", "SpectraList")
    mask_info.addColumn("double", "XMin")
    mask_info.addColumn("double", "XMax")
    for (i, value) in enumerate(qmin):
        mask_info.addRow([str(i), 0.0, value])
    for (i, value) in enumerate(qmax):
        mask_info.addRow([str(i), value, 100.0])

    return mask_info


def get_each_spectra_xmin_xmax(wksp):
    ''' Get Xmin and Xmax lists for Workspace, excluding
    values of inf and NaN

    :param wksp: Workspace to extract the xmin and xmax values from
    :type qmin: Mantid Workspace

    :return: Lists for XMin and XMax values:
    :rtype: (list, list) == (xmin, xmax)
    '''
    xmin = list()
    xmax = list()
    numSpectra = wksp.getNumberHistograms()
    for i in range(numSpectra):
        x = wksp.readX(i)
        xmin.append(np.nanmin(x[x != -np.inf]))
        y = wksp.readY(i)
        num_zeros = 0
        for ii in range(1, len(y) + 1):
            if y[-ii] != 0 and not np.isnan(y[-ii]):
                break
            num_zeros += 1
        xmax_tmp = x[len(y) - num_zeros - 1]
        xmax.append(min(xmax_tmp, np.nanmax(x[x != np.inf])))
    return xmin, xmax

# -----------------------------------------------------
# Function to expand string of ints with dashes
# Ex. "1-3, 8-9, 12" -> [1,2,3,8,9,12]


def expand_ints(s):
    spans = (el.partition('-')[::2] for el in s.split(','))
    ranges = (range(int(s), int(e) + 1 if e else int(s) + 1)
              for s, e in spans)
    all_nums = itertools.chain.from_iterable(ranges)
    return list(all_nums)

# -------------------------------------------------------------------------
# Function to compress list of ints with dashes
# Ex. [1,2,3,8,9,12] -> 1-3, 8-9, 12


def compress_ints(line_nums):
    seq = []
    final = []
    last = 0

    for index, val in enumerate(line_nums):

        if last + 1 == val or index == 0:
            seq.append(val)
            last = val
        else:
            if len(seq) > 1:
                final.append(str(seq[0]) + '-' + str(seq[len(seq) - 1]))
            else:
                final.append(str(seq[0]))
            seq = []
            seq.append(val)
            last = val

        if index == len(line_nums) - 1:
            if len(seq) > 1:
                final.append(str(seq[0]) + '-' + str(seq[len(seq) - 1]))
            else:
                final.append(str(seq[0]))

    final_str = ', '.join(map(str, final))
    return final_str

# -------------------------------------------------------------------------
# Volume in Beam


class Shape(object):
    def __init__(self):
        self.shape = None

    def getShape(self):
        return self.shape


class Cylinder(Shape):
    def __init__(self):
        self.shape = 'Cylinder'

    def volume(self, Radius=None, Height=None, **kwargs):
        return np.pi * Height * Radius * Radius


class Sphere(Shape):
    def __init__(self):
        self.shape = 'Sphere'

    def volume(self, Radius=None, **kwargs):
        return (4. / 3.) * np.pi * Radius * Radius * Radius


class GeometryFactory(object):

    @staticmethod
    def factory(Geometry):
        factory = {"Cylinder": Cylinder(),
                   "Sphere": Sphere()}
        return factory[Geometry["Shape"]]


def getNumberAtoms(PackingFraction, MassDensity, MolecularMass, Geometry=None):
    # setup the geometry of the sample
    if Geometry is None:
        Geometry = dict()
    if "Shape" not in Geometry:
        Geometry["Shape"] = 'Cylinder'

    # get sample volume in container
    space = GeometryFactory.factory(Geometry)
    volume_in_beam = space.volume(**Geometry)

    number_density = PackingFraction * MassDensity / \
        MolecularMass * Avogadro  # atoms/cm^3
    natoms = number_density * volume_in_beam  # atoms
    return natoms

# Event Filters


def GenerateEventsFilterFromFiles(filenames, OutputWorkspace,
                                  InformationWorkspace, **kwargs):

    logName = kwargs.get('LogName', None)
    minValue = kwargs.get('MinimumLogValue', None)
    maxValue = kwargs.get('MaximumLogValue', None)
    logInterval = kwargs.get('LogValueInterval', None)
    unitOfTime = kwargs.get('UnitOfTime', 'Nanoseconds')

    # TODO - handle multi-file filtering. Delete this line once implemented.
    if len(filenames) == 1:
        error = 'Multi-file filtering is not yet supported. (Stay tuned...)'
        raise Exception(error)

    for i, filename in enumerate(filenames):
        Load(Filename=filename, OutputWorkspace=filename)
        splitws, infows = GenerateEventsFilter(InputWorkspace=filename,
                                               UnitOfTime=unitOfTime,
                                               LogName=logName,
                                               MinimumLogValue=minValue,
                                               MaximumLogValue=maxValue,
                                               LogValueInterval=logInterval)
        if i == 0:
            GroupWorkspaces(splitws, OutputWorkspace=OutputWorkspace)
            GroupWorkspaces(infows, OutputWorkspace=InformationWorkspace)
        else:
            mtd[OutputWorkspace].add(splitws)
            mtd[InformationWorkspace].add(infows)
    return

# -------------------------------------------------------------------------
# Utils


def print_unit_info(workspace):
    ws = mtd[workspace]
    for i in range(ws.axes()):
        axis = ws.getAxis(i)
        print(
            "Axis {0} is a {1}{2}{3}".format(
                i,
                "Spectrum Axis" if axis.isSpectra() else "",
                "Text Axis" if axis.isText() else "",
                "Numeric Axis" if axis.isNumeric() else ""))

        unit = axis.getUnit()
        print("\n YUnit:{0}".format(ws.YUnit()))
        print("\t caption:{0}".format(unit.caption()))
        print("\t symbol:{0}".format(unit.symbol()))
    return


def SetInelasticCorrection(inelastic_dict):
    default_inelastic_dict = {"Type": None}

    if inelastic_dict is None:
        return default_inelastic_dict

    corr_type = inelastic_dict["Type"]
    if corr_type is None or corr_type == u'None':
        return default_inelastic_dict

    if corr_type:
        if corr_type == "Placzek":
            default_settings = {"Order": "1st",
                                "Self": True,
                                "Interference": True,
                                "SampleTemperature": "300",
                                "FitSpectrumWith": "GaussConvCubicSpline",
                                "LambdaBinningForFit": "0.16,0.04,2.8",
                                "LambdaBinningForCalc": "0.1,0.0001,3.0"}
            inelastic_settings = default_settings.copy()
            inelastic_settings.update(inelastic_dict)

        else:
            raise Exception("Unknown Inelastic Correction Type")

    return inelastic_settings


def get_self_scattering_level(config, max_qbinning):
    """Reads the SelfScatteringLevelCorrection option from the input config

    :param config: Input configuration dictionary
    :param max_qbinning: Maximum q binning value used to clamp the max
    level for each bank
    :return: Dictionary of bank number with tuple of min,max fit range
    or an empty dict if not specified
    """
    self_scattering_dict = dict()

    opt = "SelfScatteringLevelCorrection"
    if opt in config:
        bank_levels = config[opt]
        # return empty dict if there are no banks specified
        if len(bank_levels) == 0:
            return dict()
        for key, value in bank_levels.items():
            # get the bank number
            if not key.startswith("Bank"):
                raise RuntimeError("Expected a 'Bank' followed by number for "
                                   "SelfScatteringLevelCorrection option")
            bank = int(key.lstrip("Bank"))

            # validate the fit ranges
            if not isinstance(value, list) or len(value) != 2:
                raise RuntimeError(
                    "Expected a list of values [min, max] for each bank in "
                    "the SelfScatteringLevelCorrection option")
            if value[1] <= value[0]:
                raise RuntimeError(
                    "Max value cannot be <= min for Bank{} in "
                    "SelfScatteringLevelCorrection".format(bank))
            # clamp max to the value of Merging['QBinning'][2]
            value[1] = min(value[1], max_qbinning)
            value = tuple(value)

            self_scattering_dict[bank] = value
    return self_scattering_dict


def one_and_only_one(iterable):
    """Determine if iterable (ie list) has one and only one `True` value

    :param iterable: The iterable to check
    :type iterable: list

    :return: If there is one and only one True
    :rtype: bool
    """
    try:
        iterator = iter(iterable)
        has_true = any(iterator)
        has_another_true = any(iterator)
        return has_true and not has_another_true
    except Exception as e:
        print(e)
        raise


def find_key_match_in_dict(keys, dictionary):
    """ Check if one and only one of the keys is in the dictionary
    and return its value

    :param key: Keys we will check for in dictionary
    :type key: str
    :param dictionary: Dictionary to check
    :type dictionary: dict

    :return: Either the value in dictionary for the key or None if not found
    :rtype: value in dict or None
    """
    # Get the boolean for each key if it exists in the dictionary
    keys_exist_in_dict = map(lambda key: key in dictionary, keys)

    # If only one exists, return the match, else raise exception
    if one_and_only_one(keys_exist_in_dict):
        for key in keys:
            if key in dictionary:
                return dictionary[key]

    # None of the keys in the dictionary, return None
    return None


def extract_key_match_from_dict(keys, dictionary):
    """ Convienence function for extraction of one key from dictionary

    :param keys: Keys to check against dictionary
    :type keys: list
    :param dictionary: Dictionary to check
    "type dictionary: dict

    :return: The exctracted value
    :rtype: any
    """
    out = find_key_match_in_dict(keys, dictionary)
    if out:
        return out
    else:
        e = "No matching key found. Valid keys are {}".format(keys)
        raise Exception(e)


def get_sample(config):
    """ Extract the sample section from JSON input

    :param config: JSON input for reduction
    :type config: dict

    :return: The exctracted value for sample in the input
    :rtype: any
    """
    keys = ["Sample"]
    out = extract_key_match_from_dict(keys, config)
    return out


def get_normalization(config):
    """ Extract the normalization section from JSON input

    :param config: JSON input for reduction
    :type config: dict

    :return: The exctracted value for normalization in the input
    :rtype: any
    """
    keys = ["Normalization", "Normalisation", "Vanadium"]
    out = extract_key_match_from_dict(keys, config)
    return out


def TotalScatteringReduction(config: dict = None):
    #################################################################
    # Parse configuration from input argument 'config'
    #################################################################
    if config is None:
        raise RuntimeError('Argument config cannot be None')

    facility = config['Facility']
    title = config['Title']
    instr = config['Instrument']

    # Get an instance to Mantid's logger
    log = Logger("TotalScatteringReduction")

    # Message to be presented at the very end of the reduction via the logger.
    final_message = ''

    # Get sample info
    sample = get_sample(config)
    sam_mass_density = sample.get('MassDensity', None)
    sam_packing_fraction = sample.get('PackingFraction', None)
    sam_geometry = sample.get('Geometry', None)
    sam_material = sample.get('Material', None)

    sam_geo_dict = {'Shape': 'Cylinder',
                    'Radius': config['Sample']['Geometry']['Radius'],
                    'Height': config['Sample']['Geometry']['Height']}

    sam_eff_density = sam_mass_density * sam_packing_fraction
    sam_mat_dict = {'ChemicalFormula': sam_material,
                    'SampleMassDensity': sam_eff_density}

    if 'Environment' in config:
        sam_env_dict = {'Name': config['Environment']['Name'],
                        'Container': config['Environment']['Container']}
    else:
        sam_env_dict = {'Name': 'InAir',
                        'Container': 'PAC06'}

    abs_cache_fn = sam_mat_dict["ChemicalFormula"].replace(" ", "_").replace(".", "p")
    tmp_fn = "_md_" + "{0:7.5F}".format(sam_mat_dict['SampleMassDensity']).replace(".", "p")
    abs_cache_fn += tmp_fn
    abs_cache_fn += ("_pf_" + "{0:3.1F}".format(sam_packing_fraction).replace(".", "p"))
    abs_cache_fn += ("_sp_" + sam_geo_dict['Shape'])
    abs_cache_fn += ("_r_" + "{0:6.4F}".format(sam_geo_dict['Radius']).replace(".", "p"))
    abs_cache_fn += ("_h_" + "{0:3.1F}".format(sam_geo_dict['Height']).replace(".", "p"))
    abs_cache_fn += ("_env_" + sam_env_dict['Name'])
    abs_cache_fn += ("_cont_" + sam_env_dict['Container'])

    # Get normalization info
    van = get_normalization(config)
    van_mass_density = van.get('MassDensity', None)
    van_packing_fraction = van.get('PackingFraction', 1.0)
    van_geometry = van.get('Geometry', None)
    van_material = van.get('Material', 'V')

    van_geo_dict = {'Shape': 'Cylinder',
                    'Radius': config['Normalization']['Geometry']['Radius'],
                    'Height': config['Normalization']['Geometry']['Height']}

    van_eff_density = van_mass_density * van_packing_fraction
    van_mat_dict = {'ChemicalFormula': van_material,
                    'SampleMassDensity': van_eff_density}

    abs_cache_fn_v = van_mat_dict["ChemicalFormula"].replace(" ", "_").replace(".", "p")
    tmp_fn = "_md_" + "{0:7.5F}".format(van_mat_dict['SampleMassDensity']).replace(".", "p")
    abs_cache_fn_v += tmp_fn
    abs_cache_fn_v += ("_pf_" + "{0:3.1F}".format(van_packing_fraction).replace(".", "p"))
    abs_cache_fn_v += ("_sp_" + van_geo_dict['Shape'])
    abs_cache_fn_v += ("_r_" + "{0:6.4F}".format(van_geo_dict['Radius']).replace(".", "p"))
    abs_cache_fn_v += ("_h_" + "{0:3.1F}".format(van_geo_dict['Height']).replace(".", "p"))
    abs_cache_fn_v += ("_env_" + sam_env_dict['Name'])

    # Get calibration, characterization, and other settings
    merging = config['Merging']
    binning = merging['QBinning']
    characterizations = merging.get('Characterizations', None)

    # Get the self scattering option for each bank
    self_scattering_level_correction = get_self_scattering_level(config,
                                                                 binning[2])
    if not isinstance(self_scattering_level_correction, dict):
        raise RuntimeError()

    # Get Resonance filter configuration
    res_filter = config.get('ResonanceFilter', None)
    if res_filter is not None:
        res_filter_axis = res_filter.get('Axis', None)
        res_filter_lower = res_filter.get('LowerLimits', None)
        res_filter_upper = res_filter.get('UpperLimits', None)

    # Grouping
    grouping = merging.get('Grouping', None)
    cache_dir = config.get("CacheDir", None)
    OutputDir = config.get("OutputDir", os.path.abspath('.'))

    debug_mode = config.get("DebugMode", False)
    if debug_mode:
        print("[Info] Debug mode enabled. Intermediate workspaces to be saved.")
        print("[Info] Reduction time is supposed to increase significantly.")
    else:
        print("[Info] Debug mode disabled. Only final workspace will be saved.")

    # Vanadium peak stripping parameters
    van_ps = config.get("StripVanPeaksParams", None)
    if van_ps is not None:
        van_ps_fwhm = van_ps.get("FWHM", 7)
        van_ps_tol = van_ps.get("Tolerance", 4)
        van_ps_bkg_type = van_ps.get("BackgroundType", "Quadratic")
        van_ps_hb = van_ps.get("HighBackground", True)
        van_ps_pp_tol = van_ps.get("PeakPositionTolerance", 0.01)
    else:
        van_ps_fwhm = 7
        van_ps_tol = 4
        van_ps_bkg_type = "Quadratic"
        van_ps_hb = True
        van_ps_pp_tol = 0.01

    # Create Nexus file basenames
    sample['Runs'] = expand_ints(sample['Runs'])
    sample['Background']['Runs'] = expand_ints(
        sample['Background'].get('Runs', None))

    '''
    Currently not implemented:
    # wkspIndices = merging.get('SumBanks', None)
    # high_q_linear_fit_range = config['HighQLinearFitRange']

    POWGEN options not used
    #alignAndFocusArgs['RemovePromptPulseWidth'] = 50
    # alignAndFocusArgs['CompressTolerance'] use defaults
    # alignAndFocusArgs['UnwrapRef'] POWGEN option
    # alignAndFocusArgs['LowResRef'] POWGEN option
    # alignAndFocusArgs['LowResSpectrumOffset'] POWGEN option

    How much of each bank gets merged has info here in the form of
    # {"ID", "Qmin", "QMax"}
    # alignAndFocusArgs['CropWavelengthMin'] from characterizations file
    # alignAndFocusArgs['CropWavelengthMax'] from characterizations file
    '''

    #################################################################
    # Figure out experimental runs either with run numbers
    # and facility name or retrieve file name from 'config'
    #
    # including
    # sample, sample background,
    # container, container background,
    # vanadium, vanadium background
    #################################################################
    if facility == 'SNS':
        facility_file_format = '%s_%d'
    else:
        facility_file_format = '%s%d'

    sam_scans = ','.join([facility_file_format % (instr, num)
                          for num in sample['Runs']])
    container_scans = ','.join([facility_file_format % (instr, num)
                                for num in sample['Background']["Runs"]])
    container_bg = None
    if "Background" in sample['Background']:
        sample['Background']['Background']['Runs'] = expand_ints(
            sample['Background']['Background']['Runs'])
        container_bg = ','.join([facility_file_format % (
            instr, num) for num in sample['Background']['Background']['Runs']])
        if len(container_bg) == 0:
            container_bg = None

    van['Runs'] = expand_ints(van['Runs'])
    van_scans = ','.join([facility_file_format % (instr, num)
                          for num in van['Runs']])

    van_bg_scans = None
    if 'Background' in van:
        van_bg_scans = van['Background']['Runs']
        van_bg_scans = expand_ints(van_bg_scans)
        van_bg_scans = ','.join([facility_file_format %
                                 (instr, num) for num in van_bg_scans])

    # Override Nexus file basename with Filenames if present
    if "Filenames" in sample:
        sam_scans = ','.join(sample["Filenames"])
    if "Filenames" in sample['Background']:
        container_scans = ','.join(sample['Background']["Filenames"])
    if "Background" in sample['Background']:
        if "Filenames" in sample['Background']['Background']:
            container_bg = ','.join(
                sample['Background']['Background']['Filenames'])
    if "Filenames" in van:
        van_scans = ','.join(van["Filenames"])
    if "Background" in van:
        if "Filenames" in van['Background']:
            van_bg_scans = ','.join(van['Background']["Filenames"])

    # Output nexus filename
    nexus_filename = title + '.nxs'
    try:
        os.remove(os.path.join(OutputDir, nexus_filename))
        msg = "Old NeXusfile found. Will delete it."
        msg1 = "Old NeXus file: {}"
        log.notice(msg)
        log.notice(msg1.format(os.path.join(OutputDir, nexus_filename)))
    except OSError:
        msg = "Old NeXus file not found. Moving forward."
        log.notice(msg)
        pass

    #################################################################
    # Process absorption, multiple scattering and inelastic
    # correction setup from input 'config'
    # and
    # Create absorption workspace for
    # - sample and container
    # - vanadium
    #################################################################
    # Get sample corrections
    new_abs_methods = ['SampleOnly', 'SampleAndContainer', 'FullPaalmanPings']
    sam_abs_corr = sample.get("AbsorptionCorrection", None)
    sam_ms_corr = sample.get("MultipleScatteringCorrection", None)
    sam_inelastic_corr = SetInelasticCorrection(
        sample.get('InelasticCorrection', None))
    # get the element size
    sam_abs_ms_param = sample.get("AbsMSParameters", None)
    sam_elementsize = 1.0  # mm
    con_elementsize = 1.0  # mm
    if sam_abs_ms_param:
        elementsize = sam_abs_ms_param.get("ElementSize", 1.0)
        if type(elementsize) == list:
            sam_elementsize = elementsize[0]
            con_elementsize = elementsize[1]
        else:
            sam_elementsize = elementsize
            con_elementsize = elementsize

    # Compute the absorption correction on the sample if it was provided
    sam_abs_ws = ''
    con_abs_ws = ''
    if sam_abs_corr:
        if sam_abs_corr["Type"] in new_abs_methods:
            msg = "Applying '{}' absorption correction to sample"
            log.notice(msg.format(sam_abs_corr["Type"]))
            sam_ms_method = None
            if sam_ms_corr:
                sam_ms_method = sam_ms_corr.get("Type", None)
                if sam_ms_method is not None:
                    log.notice(
                        f"Apply {sam_ms_method} multiple scattering correction"
                        "to sample"
                    )
            ipts = GetIPTS(Instrument='NOM',
                           RunNumber=sam_scans.split(",")[0].split("_")[1])
            abs_cache_fn_s = abs_cache_fn + "_s.nxs"
            abs_cache_fn_c = abs_cache_fn + "_c.nxs"
            central_cache_f_s = os.path.join(ipts, "shared/caching",
                                             abs_cache_fn_s)
            central_cache_f_c = os.path.join(ipts, "shared/caching",
                                             abs_cache_fn_c)
            f_s_exists = os.path.exists(central_cache_f_s)
            f_c_exists = os.path.exists(central_cache_f_c)
            redo_cond_1 = not f_s_exists
            so = sam_abs_corr["Type"] == "SampleOnly"
            redo_cond_2 = f_s_exists and (not f_c_exists) and (not so)

            if redo_cond_1 or redo_cond_2:
                sam_abs_ws, con_abs_ws = create_absorption_wksp(
                    sam_scans,
                    sam_abs_corr["Type"],
                    sam_geo_dict,
                    sam_mat_dict,
                    sam_env_dict,
                    ms_method=sam_ms_method,
                    elementsize=sam_elementsize,
                    con_elementsize=con_elementsize,
                    **config)
                # save abs workspaces to cached file
                if not os.path.exists(os.path.join(ipts, "shared/caching")):
                    os.makedirs(os.path.join(ipts, "shared/caching"))
                SaveNexus(InputWorkspace=sam_abs_ws,
                          Filename=central_cache_f_s)
                if con_abs_ws != "":
                    SaveNexus(InputWorkspace=con_abs_ws,
                              Filename=central_cache_f_c)
            else:
                # read cached abs workspaces
                msg = "Cached absorption file found for sample."
                msg += " Will load and use it."
                log.notice(msg)
                sam_abs_ws = LoadNexus(Filename=central_cache_f_s)
                if f_c_exists:
                    msg = "Cached absorption file found for container."
                    msg += " Will load and use it."
                    log.notice(msg)
                    con_abs_ws = LoadNexus(Filename=central_cache_f_c)
                else:
                    con_abs_ws = ""

    # Get vanadium corrections
    van_mass_density = van.get('MassDensity', van_mass_density)
    # FIXME - van_packing_fraction is specified but not used
    van_packing_fraction = van.get(
        'PackingFraction',
        van_packing_fraction)
    van_abs_corr = van.get("AbsorptionCorrection", {"Type": None})
    van_ms_corr = van.get("MultipleScatteringCorrection", {"Type": None})
    van_inelastic_corr = SetInelasticCorrection(
        van.get('InelasticCorrection', None))
    # get the elementsize for vanadium
    van_abs_ms_param = van.get("AbsMSParameters", None)
    van_elementsize = 1.0
    if van_abs_ms_param:
        van_elementsize = van_abs_ms_param.get("ElementSize", 1.0)

    # Compute the absorption correction for the vanadium if provided
    van_abs_corr_ws = ''
    if van_abs_corr:
        if van_abs_corr["Type"] in new_abs_methods:
            msg = "Applying '{}' absorption correction to vanadium"
            log.notice(msg.format(van_abs_corr["Type"]))
            van_ms_method = None
            if van_ms_corr:
                van_ms_method = van_ms_corr.get("Type", None)
                if van_ms_method is not None:
                    log.notice(
                        f"Apply {van_ms_method} multiple scattering correction"
                        "to vanadium"
                    )
            abs_cache_fn_v += "_vanadium.nxs"
            central_cache_f_v = os.path.join("/SNS/NOM/shared/mts_abs_ms",
                                             abs_cache_fn_v)
            if not os.path.exists(central_cache_f_v):
                van_abs_corr_ws, van_con_ws = create_absorption_wksp(
                    van_scans,
                    van_abs_corr["Type"],
                    van_geo_dict,
                    van_mat_dict,
                    ms_method=van_ms_method,
                    elementsize=van_elementsize,
                    **config)
                SaveNexus(InputWorkspace=van_abs_corr_ws,
                          Filename=central_cache_f_v)
            else:
                msg = "Cached absorption file found for vanadium."
                msg += " Will load and use it."
                log.notice(msg)
                van_abs_corr_ws = LoadNexus(Filename=central_cache_f_v)

    #################################################################
    # Set up parameters for AlignAndFocus
    # and
    # Create calibration, mask and grouping workspace
    #################################################################
    alignAndFocusArgs = dict()
    alignAndFocusArgs['CalFilename'] = config['Calibration']['Filename']
    # alignAndFocusArgs['GroupFilename'] don't use
    # alignAndFocusArgs['Params'] = "0.,0.02,40."
    alignAndFocusArgs['ResampleX'] = -6000
    alignAndFocusArgs['Dspacing'] = False
    alignAndFocusArgs['PreserveEvents'] = False
    alignAndFocusArgs['MaxChunkSize'] = 0
    if cache_dir is not None:
        alignAndFocusArgs['CacheDir'] = os.path.abspath(cache_dir)
    # add resonance filter related properties
    # NOTE:
    #    the default behaivor is no filtering if not specified.
    if res_filter is not None:
        alignAndFocusArgs['ResonanceFilterUnits'] = res_filter_axis
        alignAndFocusArgs['ResonanceFilterLowerLimits'] = res_filter_lower
        alignAndFocusArgs['ResonanceFilterUpperLimits'] = res_filter_upper

    # Get any additional AlignAndFocusArgs from JSON input
    if "AlignAndFocusArgs" in config:
        otherArgs = config["AlignAndFocusArgs"]
        alignAndFocusArgs.update(otherArgs)

    # Setup grouping
    output_grouping = False
    grp_wksp = "wksp_output_group"

    if grouping:
        if 'Initial' in grouping:
            if grouping['Initial'] and not grouping['Initial'] == u'':
                alignAndFocusArgs['GroupFilename'] = grouping['Initial']
        if 'Output' in grouping:
            if grouping['Output'] and not grouping['Output'] == u'':
                output_grouping = True
                LoadDetectorsGroupingFile(InputFile=grouping['Output'],
                                          OutputWorkspace=grp_wksp)
    # If no output grouping specified, create it with Calibration Grouping
    load_grouping = False
    if not output_grouping:
        LoadDiffCal(alignAndFocusArgs['CalFilename'],
                    InstrumentName=instr,
                    WorkspaceName=grp_wksp.replace('_group', ''),
                    MakeGroupingWorkspace=True,
                    MakeCalWorkspace=False,
                    MakeMaskWorkspace=False)
        load_grouping = True

    # Setup the 6 bank method if no grouping specified
    if not grouping and not (output_grouping or load_grouping):
        CreateGroupingWorkspace(InstrumentName=instr,
                                GroupDetectorsBy='Group',
                                OutputWorkspace=grp_wksp)
        alignAndFocusArgs['GroupingWorkspace'] = grp_wksp

    #################################################################
    # Load, calibrate and diffraction focus
    # (1) sample  (2) container (3) container background
    # (4) vanadium (5) vanadium background
    #################################################################
    # TODO take out the RecalculatePCharge in the future once tested
    # Load Sample
    print("#-----------------------------------#")
    print("# Sample")
    print("#-----------------------------------#")
    sam_wksp = load(
        'sample',
        sam_scans,
        sam_geometry,
        sam_material,
        sam_mass_density,
        sam_abs_ws,
        **alignAndFocusArgs)
    sample_title = "sample_and_container"
    if debug_mode:
        save_banks(InputWorkspace=sam_wksp,
                   Filename=nexus_filename,
                   Title=sample_title,
                   OutputDir=OutputDir,
                   GroupingWorkspace=grp_wksp,
                   Binning=binning)

    sam_molecular_mass = mtd[sam_wksp].sample(
    ).getMaterial().relativeMolecularMass()
    natoms = getNumberAtoms(
        sam_packing_fraction,
        sam_mass_density,
        sam_molecular_mass,
        Geometry=sam_geometry)

    # Load Sample Container
    print("#-----------------------------------#")
    print("# Sample Container")
    print("#-----------------------------------#")
    container = load(
        'container',
        container_scans,
        absorption_wksp=con_abs_ws,
        **alignAndFocusArgs)
    if debug_mode:
        save_banks(
            InputWorkspace=container,
            Filename=nexus_filename,
            Title=container,
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)

    # -------------------------------- #
    # Load Sample Container Background #
    # -------------------------------- #
    # NOTE: sample container background IS the empty instrument
    # The full formula
    #      alpha_s(I_s - I_e) - alpha_c(I_c - I_e)
    # I = ----------------------------------------
    #          alpha_v (I_v - I_v,e)
    #
    #      alpha_s I_s - alpha_c I_c - alpha_e I_e
    #   = -----------------------------------------
    #          alpha_v I_v - alpha_v I_v,e
    # where
    #                                                 A_c - A_s
    # alpha_e = alpha_s - alpha_c = 1/A_s - 1/A_c = -------------
    #                                                 A_s * A_c
    if container_bg is not None:
        print("#-----------------------------------#")
        print("# Sample Container's Background")
        print("#-----------------------------------#")
        # NOTE: to make life easier
        # alpha_e I_e = alpha_s I_e - alpha_c I_e
        container_bg_fn = container_bg
        if sam_abs_corr is not None:
            container_bg = load(
                'container_background',
                container_bg_fn,
                absorption_wksp=sam_abs_ws,
                **alignAndFocusArgs)
            tmp = load(
                'container_background_tmp',
                container_bg_fn,
                absorption_wksp=con_abs_ws,
                **alignAndFocusArgs)
            Minus(
                LHSWorkspace=container_bg,
                RHSWorkspace=tmp,
                OutputWorkspace=container_bg)
        else:
            container_bg = load(
                'container_background',
                container_bg_fn,
                absorption_wksp=sam_abs_ws,
                **alignAndFocusArgs)
        if debug_mode:
            save_banks(
                InputWorkspace=container_bg,
                Filename=nexus_filename,
                Title=container_bg,
                OutputDir=OutputDir,
                GroupingWorkspace=grp_wksp,
                Binning=binning)

    # Load Vanadium
    print("#-----------------------------------#")
    print("# Vanadium")
    print("#-----------------------------------#")
    van_wksp = load(
        'vanadium',
        van_scans,
        van_geometry,
        van_material,
        van_mass_density,
        van_abs_corr_ws,
        **alignAndFocusArgs)

    vanadium_title = "vanadium_and_background"
    if debug_mode:
        save_banks(
            InputWorkspace=van_wksp,
            Filename=nexus_filename,
            Title=vanadium_title,
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)

    van_material = mtd[van_wksp].sample().getMaterial()
    van_molecular_mass = van_material.relativeMolecularMass()
    nvan_atoms = getNumberAtoms(
        1.0,
        van_mass_density,
        van_molecular_mass,
        Geometry=van_geometry)

    print("Sample natoms:", natoms)
    print("Vanadium natoms:", nvan_atoms)
    print("Vanadium natoms / Sample natoms:", nvan_atoms / natoms)

    # ------------------------ #
    # Load Vanadium Background #
    # ------------------------ #
    # NOTE:
    # The full formula
    #      alpha_s(I_s - I_e) - alpha_c(I_c - I_e)
    # I = ----------------------------------------
    #          alpha_v (I_v - I_v,e)
    #
    #      alpha_s I_s - alpha_c I_c - alpha_e I_e
    #   = -----------------------------------------
    #          alpha_v I_v - alpha_v I_v,e
    #
    # where
    #   * I_v,e is vanadium background
    #   * alpha_e = (alpha_s - alpha_c)
    #
    # ALSO, alpha is the inverse of [effective] absorption coefficient, i.e.
    #                alpha = 1/A
    van_bg = None
    if van_bg_scans is not None:
        print("#-----------------------------------#")
        print("# Vanadium Background")
        print("#-----------------------------------#")
        # van_bg = alpha_v I_v,e
        van_bg = load(
            'vanadium_background',
            van_bg_scans,
            absorption_wksp=van_abs_corr_ws,
            **alignAndFocusArgs)

        vanadium_bg_title = "vanadium_background"
        if debug_mode:
            save_banks(
                InputWorkspace=van_bg,
                Filename=nexus_filename,
                Title=vanadium_bg_title,
                OutputDir=OutputDir,
                GroupingWorkspace=grp_wksp,
                Binning=binning)

    # Load Instrument Characterizations
    if characterizations:
        PDDetermineCharacterizations(
            InputWorkspace=sam_wksp,
            Characterizations='characterizations',
            ReductionProperties='__snspowderreduction')
        propMan = PropertyManagerDataService.retrieve('__snspowderreduction')
        qmax = 2. * np.pi / propMan['d_min'].value
        qmin = 2. * np.pi / propMan['d_max'].value
        for a, b in zip(qmin, qmax):
            print('Qrange:', a, b)
        # TODO: Add when we apply Qmin, Qmax cropping
        # mask_info = generate_cropping_table(qmin, qmax)

    #################################################################
    # STEP 1: Subtract Backgrounds
    # van       = van - van_bg
    # sam       = sam - container
    # container = container - container_bg
    # and save (1) van (2) sam (3) container
    #################################################################
    if debug_mode:
        sam_raw = 'sam_raw'
        CloneWorkspace(
            InputWorkspace=sam_wksp,
            OutputWorkspace=sam_raw)  # for later

        container_raw = 'container_raw'
        CloneWorkspace(
            InputWorkspace=container,
            OutputWorkspace=container_raw)  # for later

    if van_bg is not None:
        RebinToWorkspace(
            WorkspaceToRebin=van_bg,
            WorkspaceToMatch=van_wksp,
            OutputWorkspace=van_bg)
        Minus(
            LHSWorkspace=van_wksp,
            RHSWorkspace=van_bg,
            OutputWorkspace=van_wksp)

    RebinToWorkspace(
        WorkspaceToRebin=container,
        WorkspaceToMatch=sam_wksp,
        OutputWorkspace=container)
    Minus(
        LHSWorkspace=sam_wksp,
        RHSWorkspace=container,
        OutputWorkspace=sam_wksp)

    # If no absorption correction is to be performed, we don't
    # need to subtract the container bkg. Refer to the note
    # in the `Load Sample Container Background` section above
    # for the mechanism being used.
    if container_bg is not None and sam_abs_corr is not None:
        RebinToWorkspace(
            WorkspaceToRebin=container_bg,
            WorkspaceToMatch=sam_wksp,
            OutputWorkspace=container_bg)
        Minus(
            LHSWorkspace=sam_wksp,
            RHSWorkspace=container_bg,
            OutputWorkspace=sam_wksp)

    container_title = "container_minus_back"
    vanadium_title = "vanadium_minus_back"
    sample_title = "sample_minus_back"
    if debug_mode:
        for wksp in [container, van_wksp, sam_wksp]:
            if mtd[wksp].getDimension(0).getUnits() != "Angstrom^-1":
                ConvertUnits(
                    InputWorkspace=wksp,
                    OutputWorkspace=wksp,
                    Target="MomentumTransfer",
                    EMode="Elastic")
        RebinToWorkspace(
            WorkspaceToRebin=container,
            WorkspaceToMatch=sam_wksp,
            OutputWorkspace=container)
        RebinToWorkspace(
            WorkspaceToRebin=container_bg,
            WorkspaceToMatch=sam_wksp,
            OutputWorkspace=container_bg)
        Minus(
            LHSWorkspace=container,
            RHSWorkspace=container_bg,
            OutputWorkspace=container)
        save_banks(
            InputWorkspace=container,
            Filename=nexus_filename,
            Title=container_title,
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)
        save_banks(
            InputWorkspace=van_wksp,
            Filename=nexus_filename,
            Title=vanadium_title,
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)
        save_banks(
            InputWorkspace=sam_wksp,
            Filename=nexus_filename,
            Title=sample_title,
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)

    #################################################################
    # STEP 2.0: Prepare vanadium as normalization calibration
    # Multiple-Scattering and Absorption (Steps 2-4) for Vanadium
    # Vanadium peak strip, smooth, Placzek
    #
    # Numerical way of absorption and multiple scattering correction
    # has been implemented at the stage of loading data.
    #
    #################################################################
    van_corrected = 'van_corrected'
    cond_v_1 = van_abs_corr['Type'] == 'Carpenter' \
        or van_ms_corr['Type'] == 'Carpenter'
    cond_v_2 = van_abs_corr['Type'] == 'Mayers' \
        or van_ms_corr['Type'] == 'Mayers'
    if cond_v_1 or cond_v_2:
        ConvertUnits(
            InputWorkspace=van_wksp,
            OutputWorkspace=van_corrected,
            Target="Wavelength",
            EMode="Elastic")

    if "Type" in van_abs_corr:
        if van_abs_corr['Type'] == 'Carpenter' \
                or van_ms_corr['Type'] == 'Carpenter':
            CarpenterSampleCorrection(
                InputWorkspace=van_corrected,
                OutputWorkspace=van_corrected,
                CylinderSampleRadius=van['Geometry']['Radius'])
        elif van_abs_corr['Type'] == 'Mayers' \
                or van_ms_corr['Type'] == 'Mayers':
            if van_ms_corr['Type'] == 'Mayers':
                MayersSampleCorrection(
                    InputWorkspace=van_corrected,
                    OutputWorkspace=van_corrected,
                    MultipleScattering=True)
            else:
                MayersSampleCorrection(
                    InputWorkspace=van_corrected,
                    OutputWorkspace=van_corrected,
                    MultipleScattering=False)
        else:
            RenameWorkspace(
                InputWorkspace=van_wksp,
                OutputWorkspace=van_corrected)
    else:
        RenameWorkspace(
            InputWorkspace=van_wksp,
            OutputWorkspace=van_corrected)

    if cond_v_1 or cond_v_2:
        ConvertUnits(
            InputWorkspace=van_corrected,
            OutputWorkspace=van_corrected,
            Target='MomentumTransfer',
            EMode='Elastic')

    vanadium_title += "_ms_abs_corrected"
    if debug_mode:
        save_banks(
            InputWorkspace=van_corrected,
            Filename=nexus_filename,
            Title=vanadium_title,
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)
        save_banks(
            InputWorkspace=van_corrected,
            Filename=nexus_filename,
            Title=vanadium_title + "_with_peaks",
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)

    # Smooth Vanadium (strip peaks plus smooth)
    ConvertUnits(
        InputWorkspace=van_corrected,
        OutputWorkspace=van_corrected,
        Target='dSpacing',
        EMode='Elastic')

    # After StripVanadiumPeaks, the workspace goes from EventWorkspace ->
    # Workspace2D
    StripVanadiumPeaks(
        InputWorkspace=van_corrected,
        OutputWorkspace=van_corrected,
        FWHM=van_ps_fwhm,
        Tolerance=van_ps_tol,
        BackgroundType=van_ps_bkg_type,
        HighBackground=van_ps_hb,
        PeakPositionTolerance=van_ps_pp_tol)

    vanadium_title += '_peaks_stripped'
    if debug_mode:
        ConvertUnits(
            InputWorkspace=van_corrected,
            OutputWorkspace=van_corrected,
            Target='MomentumTransfer',
            EMode='Elastic')
        save_banks(
            InputWorkspace=van_corrected,
            Filename=nexus_filename,
            Title=vanadium_title,
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)

    ConvertUnits(
        InputWorkspace=van_corrected,
        OutputWorkspace=van_corrected,
        Target='TOF',
        EMode='Elastic')

    FFTSmooth(
        InputWorkspace=van_corrected,
        OutputWorkspace=van_corrected,
        Filter="Butterworth",
        Params='20,2',
        IgnoreXBins=True,
        AllSpectra=True)

    ConvertUnits(
        InputWorkspace=van_corrected,
        OutputWorkspace=van_corrected,
        Target='MomentumTransfer',
        EMode='Elastic')

    if debug_mode:
        vanadium_title += '_smoothed'
        save_banks(
            InputWorkspace=van_corrected,
            Filename=nexus_filename,
            Title=vanadium_title,
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)

    # Inelastic correction
    if van_inelastic_corr['Type'] == "Placzek":
        van_scan = van['Runs'][0]
        van_incident_wksp = 'van_incident_wksp'
        van_inelastic_opts = SetInelasticCorrection(
            van.get('InelasticCorrection', None))
        lambda_binning_fit = van_inelastic_opts['LambdaBinningForFit']
        lambda_binning_calc = van_inelastic_opts['LambdaBinningForCalc']
        print('van_scan:', van_scan)

        monitor_wksp = 'monitor'
        eff_corr_wksp = 'monitor_efficiency_correction_wksp'
        LoadNexusMonitors(Filename=facility_file_format % (instr, van_scan),
                          OutputWorkspace=monitor_wksp)
        ConvertUnits(InputWorkspace=monitor_wksp, OutputWorkspace=monitor_wksp,
                     Target='Wavelength', EMode='Elastic')
        Rebin(InputWorkspace=monitor_wksp,
              OutputWorkspace=monitor_wksp,
              Params=lambda_binning_fit,
              PreserveEvents=False)
        ConvertToPointData(InputWorkspace=monitor_wksp,
                           OutputWorkspace=monitor_wksp)
        CalculateEfficiencyCorrection(WavelengthRange=lambda_binning_fit,
                                      ChemicalFormula="(He3)",
                                      DensityType="Number Density",
                                      Density=1.93138101e-08,
                                      Thickness=.1,
                                      OutputWorkspace=eff_corr_wksp)
        Multiply(
            LHSWorkspace=monitor_wksp,
            RHSWorkspace=eff_corr_wksp,
            OutputWorkspace=van_incident_wksp)
        mtd[van_incident_wksp].setDistribution(True)

        fit_type = van['InelasticCorrection']['FitSpectrumWith']
        FitIncidentSpectrum(
            InputWorkspace=van_incident_wksp,
            OutputWorkspace=van_incident_wksp,
            FitSpectrumWith=fit_type,
            BinningForFit=lambda_binning_fit,
            BinningForCalc=lambda_binning_calc,
            PlotDiagnostics=False)

        van_placzek = 'van_placzek'

        SetSample(
            InputWorkspace=van_incident_wksp,
            Material=van_mat_dict)

        calc_interfere = van['InelasticCorrection']['Interference']
        if calc_interfere:
            van_t = float(van['InelasticCorrection']['SampleTemperature'])
            van_pl_order = 2
        else:
            van_t = None
            van_pl_order = 1

        CalculatePlaczek(
            InputWorkspace=van_corrected,
            IncidentSpectra=van_incident_wksp,
            OutputWorkspace=van_placzek,
            LambdaD=1.44,
            Order=van_pl_order,
            ScaleByPackingFraction=False,
            SampleTemperature=van_t)

        ConvertToHistogram(
            InputWorkspace=van_placzek,
            OutputWorkspace=van_placzek)
        if mtd[van_corrected].YUnit() == "":
            mtd[van_corrected].setYUnit("Counts")
        mtd[van_placzek].setYUnit("Counts")

        if not mtd[van_placzek].isDistribution():
            ConvertToDistribution(van_placzek)

        bin_tmp = float(lambda_binning_calc.split(',')[1])
        mtd[van_placzek] = mtd[van_placzek] * bin_tmp

        # Rebin and save in Q
        ConvertUnits(
            InputWorkspace=van_placzek,
            OutputWorkspace=van_placzek,
            Target='MomentumTransfer',
            EMode='Elastic')

        Rebin(
            InputWorkspace=van_placzek,
            OutputWorkspace=van_placzek,
            Params=binning,
            PreserveEvents=True)

        if debug_mode:
            save_banks(
                InputWorkspace=van_placzek,
                Filename=nexus_filename,
                Title="vanadium_placzek",
                OutputDir=OutputDir,
                GroupingWorkspace=grp_wksp,
                Binning=binning)

    SetUncertainties(
        InputWorkspace=van_corrected,
        OutputWorkspace=van_corrected,
        SetError='zero')

    #################################################################
    # STEP 2.1: Normalize by Vanadium and
    #           convert to Q (momentum transfer)
    #           For sample, container, sample raw, vanadium background
    #################################################################
    if debug_mode:
        wksp_list = [sam_wksp, sam_raw, van_corrected]
    else:
        wksp_list = [sam_wksp, van_corrected]

    for name in wksp_list:
        if not debug_mode and name is not van_corrected:
            if mtd[name].getDimension(0).getUnits() != "Angstrom^-1":
                ConvertUnits(
                    InputWorkspace=name,
                    OutputWorkspace=name,
                    Target='MomentumTransfer',
                    EMode='Elastic',
                    ConvertFromPointData=False)

        Rebin(
            InputWorkspace=name,
            OutputWorkspace=name,
            Params=binning,
            PreserveEvents=True)

    # Save (sample - back) / van_corrected
    Divide(
        LHSWorkspace=sam_wksp,
        RHSWorkspace=van_corrected,
        OutputWorkspace=sam_wksp)

    sample_title += "_normalized"
    if debug_mode:
        save_banks(
            InputWorkspace=sam_wksp,
            Filename=nexus_filename,
            Title=sample_title,
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)

        # Save the sample / normalized (i.e., no background subtraction)
        Divide(
           LHSWorkspace=sam_raw,
           RHSWorkspace=van_corrected,
           OutputWorkspace=sam_raw)

        save_banks(
            InputWorkspace=sam_raw,
            Filename=nexus_filename,
            Title="sample_normalized",
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)

    # Output an initial I(Q) for sample
    if debug_mode:
        iq_filename = title + '_initial_iofq_banks.nxs'
        try:
            os.remove(os.path.join(OutputDir, iq_filename))
            msg = "Old NeXus file found for initial iofq. Will delete it."
            msg1 = "Old NeXus file: {}"
            log.notice(msg)
            log.notice(msg1.format(os.path.join(OutputDir, iq_filename)))
        except OSError:
            msg = "Old NeXus file not found for initial iofq. Moving forward."
            log.notice(msg)
            pass
        save_banks(
            InputWorkspace=sam_wksp,
            Filename=iq_filename,
            Title="IQ_banks",
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)

        wksp_list = [container, container_raw, van_corrected]

        if container_bg is not None:
            wksp_list.append(container_bg)
        if van_bg is not None:
            wksp_list.append(van_bg)

        for name in wksp_list:
            if mtd[name].getDimension(0).getUnits() != "Angstrom^-1":
                ConvertUnits(
                    InputWorkspace=name,
                    OutputWorkspace=name,
                    Target='MomentumTransfer',
                    EMode='Elastic',
                    ConvertFromPointData=False)

            Rebin(
                InputWorkspace=name,
                OutputWorkspace=name,
                Params=binning,
                PreserveEvents=True)

        # Save the container - container_background / normalized
        Divide(
            LHSWorkspace=container,
            RHSWorkspace=van_corrected,
            OutputWorkspace=container)

        container_title += '_normalized'
        save_banks(
            InputWorkspace=container,
            Filename=nexus_filename,
            Title=container_title,
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)

        # Save the container / normalized (ie no background subtraction)
        Divide(
           LHSWorkspace=container_raw,
           RHSWorkspace=van_corrected,
           OutputWorkspace=container_raw)

        save_banks(
            InputWorkspace=container_raw,
            Filename=nexus_filename,
            Title="container_normalized",
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)

        # Save the container_background / normalized
        if container_bg is not None:
            Divide(
                LHSWorkspace=container_bg,
                RHSWorkspace=van_corrected,
                OutputWorkspace=container_bg)

            container_bg_title = "container_back_normalized"
            save_banks(
                InputWorkspace=container_bg,
                Filename=nexus_filename,
                Title=container_bg_title,
                OutputDir=OutputDir,
                GroupingWorkspace=grp_wksp,
                Binning=binning)

        # Save the vanadium_background / normalized
        if van_bg is not None:
            Divide(
                LHSWorkspace=van_bg,
                RHSWorkspace=van_corrected,
                OutputWorkspace=van_bg)

            vanadium_bg_title += "_normalized"
            save_banks(
                InputWorkspace=van_bg,
                Filename=nexus_filename,
                Title=vanadium_bg_title,
                OutputDir=OutputDir,
                GroupingWorkspace=grp_wksp,
                Binning=binning)

    #################################################################
    # STEP 3 & 4: Subtract multiple scattering and apply absorption correction
    #################################################################
    sam_corrected = 'sam_corrected'
    if sam_abs_corr and sam_ms_corr:
        # When full PP approach is used for absorption correction, the
        # calculation was performed when loading data at the stage of
        # align and focus. Multiple scattering calculation was also
        # embedded there.
        cond_1 = sam_abs_corr['Type'] == 'Carpenter' \
            or sam_ms_corr['Type'] == 'Carpenter'
        cond_2 = sam_abs_corr['Type'] == 'Mayers' \
            or sam_ms_corr['Type'] == 'Mayers'
        if cond_1 or cond_2:
            ConvertUnits(
                InputWorkspace=sam_wksp,
                OutputWorkspace=sam_wksp,
                Target="Wavelength",
                EMode="Elastic")

        if sam_abs_corr['Type'] == 'Carpenter' \
                or sam_ms_corr['Type'] == 'Carpenter':
            CarpenterSampleCorrection(
                InputWorkspace=sam_wksp,
                OutputWorkspace=sam_corrected,
                CylinderSampleRadius=sample['Geometry']['Radius'])
        elif sam_abs_corr['Type'] == 'Mayers' \
                or sam_ms_corr['Type'] == 'Mayers':
            if sam_ms_corr['Type'] == 'Mayers':
                MayersSampleCorrection(
                    InputWorkspace=sam_wksp,
                    OutputWorkspace=sam_corrected,
                    MultipleScattering=True)
            else:
                MayersSampleCorrection(
                    InputWorkspace=sam_wksp,
                    OutputWorkspace=sam_corrected,
                    MultipleScattering=False)
        else:
            RenameWorkspace(
                InputWorkspace=sam_wksp,
                OutputWorkspace=sam_corrected)

        if cond_1 or cond_2:
            ConvertUnits(
                InputWorkspace=sam_corrected,
                OutputWorkspace=sam_corrected,
                Target='MomentumTransfer',
                EMode='Elastic')

        sample_title += "_ms_abs_corrected"
        if debug_mode:
            save_banks(
                InputWorkspace=sam_corrected,
                Filename=nexus_filename,
                Title=sample_title,
                OutputDir=OutputDir,
                GroupingWorkspace=grp_wksp,
                Binning=binning)
    else:
        RenameWorkspace(InputWorkspace=sam_wksp,
                        OutputWorkspace=sam_corrected)

    #################################################################
    # STEP 5: Divide by number of atoms in sample
    #################################################################

    mtd[sam_corrected] = (nvan_atoms / natoms) * mtd[sam_corrected]
    ConvertUnits(
        InputWorkspace=sam_corrected,
        OutputWorkspace=sam_corrected,
        Target='MomentumTransfer',
        EMode='Elastic')

    sample_title += "_norm_by_atoms"
    if debug_mode:
        save_banks(
            InputWorkspace=sam_corrected,
            Filename=nexus_filename,
            Title=sample_title,
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)

    #################################################################
    # STEP 6: Divide by total scattering length squared = total scattering
    # cross-section over 4 * pi
    #################################################################
    van_material = mtd[van_corrected].sample().getMaterial()
    sigma_v = van_material.totalScatterXSection()
    prefactor = (sigma_v / (4. * np.pi))
    msg = "Total scattering cross-section of Vanadium:{} sigma_v / 4*pi: {}"
    print(msg.format(sigma_v, prefactor))

    if van_inelastic_corr['Type'] == "Placzek":
        mtd[sam_corrected] = (prefactor + mtd[van_placzek]) * mtd[sam_corrected]
    else:
        mtd[sam_corrected] = prefactor * mtd[sam_corrected]
    sample_title += '_multiply_by_vanSelfScat'
    if debug_mode:
        save_banks(
            InputWorkspace=sam_corrected,
            Filename=nexus_filename,
            Title=sample_title,
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)

    #################################################################
    # STEP 7: Inelastic correction
    #################################################################
    if sam_inelastic_corr['Type'] == "Placzek":
        if sam_material is None:
            error = "For Placzek correction, must specifiy a sample material."
            raise Exception(error)

        sam_scan = sample['Runs'][0]
        sam_incident_wksp = 'sam_incident_wksp'
        sam_inelastic_opts = SetInelasticCorrection(
            sample.get('InelasticCorrection', None))
        lambda_binning_fit = sam_inelastic_opts['LambdaBinningForFit']
        lambda_binning_calc = sam_inelastic_opts['LambdaBinningForCalc']

        monitor_wksp = 'monitor'
        eff_corr_wksp = 'monitor_efficiency_correction_wksp'
        LoadNexusMonitors(Filename=facility_file_format % (instr, sam_scan),
                          OutputWorkspace=monitor_wksp)
        ConvertUnits(InputWorkspace=monitor_wksp, OutputWorkspace=monitor_wksp,
                     Target='Wavelength', EMode='Elastic')
        Rebin(InputWorkspace=monitor_wksp,
              OutputWorkspace=monitor_wksp,
              Params=lambda_binning_fit,
              PreserveEvents=False)
        ConvertToPointData(InputWorkspace=monitor_wksp,
                           OutputWorkspace=monitor_wksp)
        CalculateEfficiencyCorrection(WavelengthRange=lambda_binning_fit,
                                      ChemicalFormula="(He3)",
                                      DensityType="Number Density",
                                      Density=1.93138101e-08,
                                      Thickness=.1,
                                      OutputWorkspace=eff_corr_wksp)
        Multiply(
            LHSWorkspace=monitor_wksp,
            RHSWorkspace=eff_corr_wksp,
            OutputWorkspace=sam_incident_wksp)
        mtd[sam_incident_wksp].setDistribution(True)

        fit_type = sample['InelasticCorrection']['FitSpectrumWith']
        FitIncidentSpectrum(
            InputWorkspace=sam_incident_wksp,
            OutputWorkspace=sam_incident_wksp,
            FitSpectrumWith=fit_type,
            BinningForFit=lambda_binning_fit,
            BinningForCalc=lambda_binning_calc)

        sam_placzek = 'sam_placzek'
        SetSample(
            InputWorkspace=sam_incident_wksp,
            Material=sam_mat_dict)

        calc_interfere = sample['InelasticCorrection']['Interference']
        if calc_interfere:
            sample_t = float(sample['InelasticCorrection']['SampleTemperature'])
            sample_pl_order = 2
        else:
            sample_t = None
            sample_pl_order = 1

        CalculatePlaczek(
            InputWorkspace=sam_corrected,
            IncidentSpectra=sam_incident_wksp,
            OutputWorkspace=sam_placzek,
            LambdaD=1.44,
            Order=sample_pl_order,
            ScaleByPackingFraction=False,
            SampleTemperature=sample_t)

        ConvertToHistogram(
            InputWorkspace=sam_placzek,
            OutputWorkspace=sam_placzek)
        if mtd[sam_corrected].YUnit() == "":
            mtd[sam_corrected].setYUnit("Counts")
        mtd[sam_placzek].setYUnit("Counts")

        if not mtd[sam_placzek].isDistribution():
            ConvertToDistribution(sam_placzek)

        bin_tmp = float(lambda_binning_calc.split(',')[1])
        mtd[sam_placzek] = mtd[sam_placzek] * bin_tmp

        ConvertUnits(
            InputWorkspace=sam_placzek,
            OutputWorkspace=sam_placzek,
            Target='MomentumTransfer',
            EMode='Elastic')

        Rebin(
            InputWorkspace=sam_placzek,
            OutputWorkspace=sam_placzek,
            Params=binning,
            PreserveEvents=True)

        if debug_mode:
            save_banks(
                InputWorkspace=sam_placzek,
                Filename=nexus_filename,
                Title="sample_placzek",
                OutputDir=OutputDir,
                GroupingWorkspace=grp_wksp,
                Binning=binning)

        Minus(
            LHSWorkspace=sam_corrected,
            RHSWorkspace=sam_placzek,
            OutputWorkspace=sam_corrected)

        sample_title += '_placzek_corrected'
        if debug_mode:
            save_banks(
                InputWorkspace=sam_corrected,
                Filename=nexus_filename,
                Title=sample_title,
                OutputDir=OutputDir,
                GroupingWorkspace=grp_wksp,
                Binning=binning)

    # TODO Since we already went from Event -> 2D workspace, can't use this
    # anymore
    print('sam:', mtd[sam_corrected].id())
    print('van:', mtd[van_corrected].id())
    if alignAndFocusArgs['PreserveEvents']:
        CompressEvents(
            InputWorkspace=sam_corrected,
            OutputWorkspace=sam_corrected)

    #################################################################
    # STEP 8:  S(Q) and F(Q), bank-by-bank
    # processing includes
    #   1. to absolute scale
    #   2. save to S(Q)
    #   3. self scattering correction and save S(Q) corrected
    #   4. save to F(Q)
    #################################################################
    if self_scattering_level_correction:
        offset, slope = \
            calculate_and_apply_fitted_levels(sam_corrected,
                                              self_scattering_level_correction)
    else:
        offset = None

    sam_corrected_norm = sam_corrected + '_norm'
    to_absolute_scale(sam_corrected, sam_corrected_norm)
    save_banks(
        InputWorkspace=sam_corrected_norm,
        Filename=nexus_filename,
        Title="SQ_banks_normalized",
        OutputDir=OutputDir,
        GroupingWorkspace=grp_wksp,
        Binning=binning)

    fq_banks = 'FQ_banks'
    to_f_of_q(sam_corrected_norm, fq_banks)
    save_banks(
        InputWorkspace=fq_banks,
        Filename=nexus_filename,
        Title="FQ_banks",
        OutputDir=OutputDir,
        GroupingWorkspace=grp_wksp,
        Binning=binning)

    # Print log information
    material = Material(sam_corrected)
    log_file_out = open(os.path.join(os.path.abspath(OutputDir),
                                     f'{title}.log'),
                        "w")
    sep_line = "==============================="
    sep_line += "=================================\n"
    log_file_out.write(sep_line)
    log_file_out.write("{0:32s}{1:<20.10F}\n".format('<b>^2:',
                                                     material.bcoh_avg_sqrd))
    log_file_out.write("{0:32s}{1:<20.10F}\n".format("<b^2>:",
                                                     material.btot_sqrd_avg))
    laue_term = material.laue_monotonic_diffuse_scat
    log_file_out.write("{0:32s}{1:<20.10F}\n".format("Laue term:", laue_term))
    tsxsection = material.totalScatterXSection()
    log_file_out.write("{0:32s}{1:<20.10F}\n".format("Sample total xsection:",
                                                     tsxsection))
    tsxsection = Material(van_corrected).totalScatterXSection()
    log_file_out.write("{0:32s}{1:<20.10F}\n".format("Vanadium total xsection:",
                                                     tsxsection))
    log_file_out.write(sep_line)

    if offset:
        header_line = "Bank    HighQ offset    HighQ Slope    Input PF"
        header_line += "    Suggested PF\n"
        log_file_out.write(header_line)
        i_pf = sam_packing_fraction
        s_pf_all = dict()
        for key, item in offset.items():
            log_file_out.write("{0:<8d}{1:<16.3F}{2:<15.3F}".format(key,
                                                                    item,
                                                                    slope[key]))
            s_pf = sam_packing_fraction * item / material.btot_sqrd_avg
            s_pf_all[key] = s_pf
            log_file_out.write("{0:<12.1F}{1:<12.1F}".format(i_pf, s_pf))
            log_file_out.write("\n")
        sep_line1 = "-------------------------------"
        sep_line1 += "---------------------------------\n"
        log_file_out.write(sep_line1)
        log_file_out.write("{0:39s}".format("Effective Val:"))
        if len(s_pf_all) == 1:
            used_s_pf = s_pf_all[1]
        else:
            used_s_pf = s_pf_all[5]
        log_file_out.write("{0:<12.1F}{1:<12.1F}\n".format(i_pf, used_s_pf))
        log_file_out.write(sep_line)

    log_file_out.close()

    #################################################################
    #   STEP 8.1 Output Bragg Diffraction
    #            convert to TOF and save to GSAS
    #################################################################
    ConvertUnits(
        InputWorkspace=sam_corrected,
        OutputWorkspace=sam_corrected,
        Target="TOF",
        EMode="Elastic")

    ConvertToHistogram(
        InputWorkspace=sam_corrected,
        OutputWorkspace=sam_corrected)

    xmin, xmax = get_each_spectra_xmin_xmax(mtd[sam_corrected])

    CropWorkspaceRagged(
        InputWorkspace=sam_corrected,
        OutputWorkspace=sam_corrected,
        Xmin=xmin,
        Xmax=xmax)

    xmin_rebin = min(xmin)
    if "TMin" in alignAndFocusArgs.keys():
        tmin = alignAndFocusArgs["TMin"]
        info_part1 = f"[Info] 'TMin = {tmin}' found in the input config file."
        info_part2 = " Will use it for Bragg output."
        print(info_part1 + info_part2)
        xmin_rebin = tmin
    xmax_rebin = max(xmax)
    # Note: For the moment, bin size for Bragg output is hard coded.
    # May need to make it user input if necessary.
    tof_binning = "{xmin},-0.0008,{xmax}".format(xmin=xmin_rebin,
                                                 xmax=xmax_rebin)

    Rebin(
        InputWorkspace=sam_corrected,
        OutputWorkspace=sam_corrected,
        Params=tof_binning)

    SaveGSS(
        InputWorkspace=sam_corrected,
        Filename=os.path.join(os.path.abspath(OutputDir), 'GSAS', title+".gsa"),
        SplitFiles=False,
        Append=False,
        MultiplyByBinWidth=True,
        Format="SLOG",
        ExtendedHeader=True)

    SaveFocusedXYE(
        InputWorkspace=sam_corrected,
        Filename=os.path.join(os.path.abspath(OutputDir),
                              'Topas',
                              title+".xye"))

    if final_message:
        log.warning(final_message)

    return mtd[sam_corrected]
