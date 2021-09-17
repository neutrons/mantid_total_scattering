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
    SetSample, \
    SetUncertainties, \
    StripVanadiumPeaks

from total_scattering.file_handling.load import load, create_absorption_wksp
from total_scattering.file_handling.save import save_banks
from total_scattering.inelastic.placzek import \
    CalculatePlaczekSelfScattering, \
    FitIncidentSpectrum, \
    GetIncidentSpectrumFromMonitor
from total_scattering.reduction.normalizations import (
    Material, calculate_and_apply_fitted_levels, to_absolute_scale, to_f_of_q)


# Constants
MS_AND_ABS_CORR_WARNING = (
    "The multiple scattering option should not be used with the "
    "absorption correction set")


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
        xmax.append(np.nanmax(x[x != np.inf]))
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
                                "LambdaBinning": "0.16,0.04,2.8"}
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


def TotalScatteringReduction(config=None):
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
    sam_mat_dict = {'ChemicalFormula': sam_material,
                    'SampleMassDensity': sam_mass_density}
    if 'Environment' in config:
        sam_env_dict = {'Name': config['Environment']['Name'],
                        'Container': config['Environment']['Container']}
    else:
        sam_env_dict = {'Name': 'InAir',
                        'Container': 'PAC06'}
    # Get normalization info
    van = get_normalization(config)
    van_mass_density = van.get('MassDensity', None)
    van_packing_fraction = van.get('PackingFraction', 1.0)
    van_geometry = van.get('Geometry', None)
    van_material = van.get('Material', 'V')

    van_geo_dict = {'Shape': 'Cylinder',
                    'Radius': config['Normalization']['Geometry']['Radius'],
                    'Height': config['Normalization']['Geometry']['Height']}
    van_mat_dict = {'ChemicalFormula': van_material,
                    'SampleMassDensity': van_mass_density}

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
    cache_dir = config.get("CacheDir", os.path.abspath('.'))
    OutputDir = config.get("OutputDir", os.path.abspath('.'))

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
        os.remove(nexus_filename)
    except OSError:
        pass

    # Get sample corrections
    sam_abs_corr = sample.get("AbsorptionCorrection", None)
    sam_ms_corr = sample.get("MultipleScatteringCorrection", None)
    sam_inelastic_corr = SetInelasticCorrection(
        sample.get('InelasticCorrection', None))

    # Warn about having absorption correction and multiple scat correction set
    if sam_abs_corr and sam_ms_corr:
        log.warning(MS_AND_ABS_CORR_WARNING)

    # Compute the absorption correction on the sample if it was provided
    sam_abs_ws = ''
    con_abs_ws = ''
    if sam_abs_corr:
        msg = "Applying '{}' absorption correction to sample"
        log.notice(msg.format(sam_abs_corr["Type"]))
        sam_abs_ws, con_abs_ws = create_absorption_wksp(
            sam_scans,
            sam_abs_corr["Type"],
            sam_geo_dict,
            sam_mat_dict,
            sam_env_dict,
            **config)

    # Get vanadium corrections
    van_mass_density = van.get('MassDensity', van_mass_density)
    van_packing_fraction = van.get(
        'PackingFraction',
        van_packing_fraction)
    van_abs_corr = van.get("AbsorptionCorrection", {"Type": None})
    van_ms_corr = van.get("MultipleScatteringCorrection", {"Type": None})
    van_inelastic_corr = SetInelasticCorrection(
        van.get('InelasticCorrection', None))

    # Warn about having absorption correction and multiple scat correction set
    if van_abs_corr["Type"] and van_ms_corr["Type"]:
        log.warning(MS_AND_ABS_CORR_WARNING)

    # Compute the absorption correction for the vanadium if provided
    van_abs_corr_ws = ''
    if van_abs_corr:
        msg = "Applying '{}' absorption correction to vanadium"
        log.notice(msg.format(van_abs_corr["Type"]))
        van_abs_corr_ws, van_con_ws = create_absorption_wksp(
            van_scans,
            van_abs_corr["Type"],
            van_geo_dict,
            van_mat_dict,
            **config)

    alignAndFocusArgs = dict()
    alignAndFocusArgs['CalFilename'] = config['Calibration']['Filename']
    # alignAndFocusArgs['GroupFilename'] don't use
    # alignAndFocusArgs['Params'] = "0.,0.02,40."
    alignAndFocusArgs['ResampleX'] = -6000
    alignAndFocusArgs['Dspacing'] = False
    alignAndFocusArgs['PreserveEvents'] = False
    alignAndFocusArgs['MaxChunkSize'] = 8
    alignAndFocusArgs['CacheDir'] = os.path.abspath(cache_dir)
    # add resonance filter related properties
    # NOTE:
    #    the default behaivor is no filtering if not specified.
    if res_filter is not None:
        alignAndFocusArgs['ResonanceFilterAxis'] = res_filter_axis
        alignAndFocusArgs['LowerLimits'] = res_filter_lower
        alignAndFocusArgs['UpperLimits'] = res_filter_upper

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
    if not output_grouping:
        LoadDiffCal(alignAndFocusArgs['CalFilename'],
                    InstrumentName=instr,
                    WorkspaceName=grp_wksp.replace('_group', ''),
                    MakeGroupingWorkspace=True,
                    MakeCalWorkspace=False,
                    MakeMaskWorkspace=False)

    # Setup the 6 bank method if no grouping specified
    if not grouping:
        CreateGroupingWorkspace(InstrumentName=instr,
                                GroupDetectorsBy='Group',
                                OutputWorkspace=grp_wksp)
        alignAndFocusArgs['GroupingWorkspace'] = grp_wksp

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
    save_banks(
        InputWorkspace=container,
        Filename=nexus_filename,
        Title=container,
        OutputDir=OutputDir,
        GroupingWorkspace=grp_wksp,
        Binning=binning)

    # Load Sample Container Background

    if container_bg is not None:
        print("#-----------------------------------#")
        print("# Sample Container's Background")
        print("#-----------------------------------#")
        container_bg = load(
            'container_background',
            container_bg,
            **alignAndFocusArgs)
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

    # Load Vanadium Background
    van_bg = None
    if van_bg_scans is not None:
        print("#-----------------------------------#")
        print("# Vanadium Background")
        print("#-----------------------------------#")
        van_bg = load('vanadium_background', van_bg_scans, **alignAndFocusArgs)
        vanadium_bg_title = "vanadium_background"
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

    # STEP 1: Subtract Backgrounds

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

    if container_bg is not None:
        RebinToWorkspace(
            WorkspaceToRebin=container_bg,
            WorkspaceToMatch=container,
            OutputWorkspace=container_bg)
        Minus(
            LHSWorkspace=container,
            RHSWorkspace=container_bg,
            OutputWorkspace=container)

    for wksp in [container, van_wksp, sam_wksp]:
        ConvertUnits(
            InputWorkspace=wksp,
            OutputWorkspace=wksp,
            Target="MomentumTransfer",
            EMode="Elastic")
    container_title = "container_minus_back"
    vanadium_title = "vanadium_minus_back"
    sample_title = "sample_minus_back"
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

    # STEP 2.0: Prepare vanadium as normalization calibrant

    # Multiple-Scattering and Absorption (Steps 2-4) for Vanadium

    van_corrected = 'van_corrected'
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
            print("NO VANADIUM absorption or multiple scattering!")
    else:
        CloneWorkspace(
            InputWorkspace=van_corrected,
            OutputWorkspace=van_corrected)

    ConvertUnits(
        InputWorkspace=van_corrected,
        OutputWorkspace=van_corrected,
        Target='MomentumTransfer',
        EMode='Elastic')
    vanadium_title += "_ms_abs_corrected"
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

    # TODO subtract self-scattering of vanadium (According to Eq. 7 of Howe,
    # McGreevey, and Howells, JPCM, 1989)

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
        BackgroundType='Quadratic')
    ConvertUnits(
        InputWorkspace=van_corrected,
        OutputWorkspace=van_corrected,
        Target='MomentumTransfer',
        EMode='Elastic')
    vanadium_title += '_peaks_stripped'
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
        van_inelastic_opts = van['InelasticCorrection']
        lambda_binning_fit = van_inelastic_opts['LambdaBinningForFit']
        lambda_binning_calc = van_inelastic_opts['LambdaBinningForCalc']
        print('van_scan:', van_scan)
        GetIncidentSpectrumFromMonitor(
            Filename=facility_file_format % (instr, van_scan),
            OutputWorkspace=van_incident_wksp)

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
            Material={'ChemicalFormula': str(van_material),
                      'SampleMassDensity': str(van_mass_density)})

        calc_interfere = van['InelasticCorrection']['Interference']
        sample_t = float(van['InelasticCorrection']['SampleTemperature'])
        CalculatePlaczekSelfScattering(
            IncidentWorkspace=van_incident_wksp,
            ParentWorkspace=van_corrected,
            OutputWorkspace=van_placzek,
            L1=19.5,
            L2=alignAndFocusArgs['L2'],
            Polar=alignAndFocusArgs['Polar'],
            CalcInterfere=calc_interfere,
            SampleT=sample_t)

        ConvertToHistogram(
            InputWorkspace=van_placzek,
            OutputWorkspace=van_placzek)

        # Save before rebin in Q
        for wksp in [van_placzek, van_corrected]:
            ConvertUnits(
                InputWorkspace=wksp,
                OutputWorkspace=wksp,
                Target='MomentumTransfer',
                EMode='Elastic')

            Rebin(
                InputWorkspace=wksp,
                OutputWorkspace=wksp,
                Params=binning,
                PreserveEvents=True)

        save_banks(
            InputWorkspace=van_placzek,
            Filename=nexus_filename,
            Title="vanadium_placzek",
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)

        # Rebin in Wavelength
        for wksp in [van_placzek, van_corrected]:
            ConvertUnits(
                InputWorkspace=wksp,
                OutputWorkspace=wksp,
                Target='Wavelength',
                EMode='Elastic')
            Rebin(
                InputWorkspace=wksp,
                OutputWorkspace=wksp,
                Params=lambda_binning_calc,
                PreserveEvents=True)

        # Save after rebin in Q
        for wksp in [van_placzek, van_corrected]:
            ConvertUnits(
                InputWorkspace=wksp,
                OutputWorkspace=wksp,
                Target='MomentumTransfer',
                EMode='Elastic')

        # Subtract correction in Wavelength
        for wksp in [van_placzek, van_corrected]:
            ConvertUnits(
                InputWorkspace=wksp,
                OutputWorkspace=wksp,
                Target='Wavelength',
                EMode='Elastic')
            if not mtd[wksp].isDistribution():
                ConvertToDistribution(wksp)

        Minus(
            LHSWorkspace=van_corrected,
            RHSWorkspace=van_placzek,
            OutputWorkspace=van_corrected)

        # Save after subtraction
        for wksp in [van_placzek, van_corrected]:
            ConvertUnits(
                InputWorkspace=wksp,
                OutputWorkspace=wksp,
                Target='MomentumTransfer',
                EMode='Elastic')

        vanadium_title += '_placzek_corrected'
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
        Target='MomentumTransfer',
        EMode='Elastic')

    SetUncertainties(
        InputWorkspace=van_corrected,
        OutputWorkspace=van_corrected,
        SetError='zero')

    # STEP 2.1: Normalize by Vanadium

    wksp_list = [sam_wksp, sam_raw, van_corrected]
    for name in wksp_list:
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

    # Save the sample - back / normalized
    Divide(
        LHSWorkspace=sam_wksp,
        RHSWorkspace=van_corrected,
        OutputWorkspace=sam_wksp)

    sample_title += "_normalized"
    save_banks(
        InputWorkspace=sam_wksp,
        Filename=nexus_filename,
        Title=sample_title,
        OutputDir=OutputDir,
        GroupingWorkspace=grp_wksp,
        Binning=binning)

    # Save the sample / normalized (ie no background subtraction)
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
    iq_filename = title + '_initial_iofq_banks.nxs'
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

    # STEP 3 & 4: Subtract multiple scattering and apply absorption correction

    ConvertUnits(
        InputWorkspace=sam_wksp,
        OutputWorkspace=sam_wksp,
        Target="Wavelength",
        EMode="Elastic")

    sam_corrected = 'sam_corrected'
    if sam_abs_corr and sam_ms_corr:
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
            print("NO SAMPLE absorption or multiple scattering!")
            CloneWorkspace(
                InputWorkspace=sam_wksp,
                OutputWorkspace=sam_corrected)

        ConvertUnits(
            InputWorkspace=sam_corrected,
            OutputWorkspace=sam_corrected,
            Target='MomentumTransfer',
            EMode='Elastic')

        sample_title += "_ms_abs_corrected"
        save_banks(
            InputWorkspace=sam_corrected,
            Filename=nexus_filename,
            Title=sample_title,
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)
    else:
        CloneWorkspace(InputWorkspace=sam_wksp, OutputWorkspace=sam_corrected)

    # STEP 5: Divide by number of atoms in sample

    mtd[sam_corrected] = (nvan_atoms / natoms) * mtd[sam_corrected]
    ConvertUnits(
        InputWorkspace=sam_corrected,
        OutputWorkspace=sam_corrected,
        Target='MomentumTransfer',
        EMode='Elastic')

    sample_title += "_norm_by_atoms"
    save_banks(
        InputWorkspace=sam_corrected,
        Filename=nexus_filename,
        Title=sample_title,
        OutputDir=OutputDir,
        GroupingWorkspace=grp_wksp,
        Binning=binning)

    # STEP 6: Divide by total scattering length squared = total scattering
    # cross-section over 4 * pi
    van_material = mtd[van_corrected].sample().getMaterial()
    sigma_v = van_material.totalScatterXSection()
    prefactor = (sigma_v / (4. * np.pi))
    msg = "Total scattering cross-section of Vanadium:{} sigma_v / 4*pi: {}"
    print(msg.format(sigma_v, prefactor))

    mtd[sam_corrected] = prefactor * mtd[sam_corrected]
    sample_title += '_multiply_by_vanSelfScat'
    save_banks(
        InputWorkspace=sam_corrected,
        Filename=nexus_filename,
        Title=sample_title,
        OutputDir=OutputDir,
        GroupingWorkspace=grp_wksp,
        Binning=binning)

    # STEP 7: Inelastic correction

    if sam_inelastic_corr['Type'] == "Placzek":
        if sam_material is None:
            error = "For Placzek correction, must specifiy a sample material."
            raise Exception(error)

        ConvertUnits(
            InputWorkspace=sam_corrected,
            OutputWorkspace=sam_corrected,
            Target='Wavelength',
            EMode='Elastic')

        for sam_scan in sample['Runs']:
            sam_incident_wksp = 'sam_incident_wksp'
            sam_inelastic_opts = sample['InelasticCorrection']
            lambda_binning_fit = sam_inelastic_opts['LambdaBinningForFit']
            lambda_binning_calc = sam_inelastic_opts['LambdaBinningForCalc']
            GetIncidentSpectrumFromMonitor(
                Filename=facility_file_format % (instr, sam_scan),
                OutputWorkspace=sam_incident_wksp)

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
                Material={'ChemicalFormula': str(sam_material),
                          'SampleMassDensity': str(sam_mass_density)})

            calc_interfere = sample['InelasticCorrection']['Interference']
            sample_t = float(sample['InelasticCorrection']['SampleTemperature'])
            CalculatePlaczekSelfScattering(
                IncidentWorkspace=sam_incident_wksp,
                ParentWorkspace=sam_corrected,
                OutputWorkspace=sam_placzek,
                L1=19.5,
                L2=alignAndFocusArgs['L2'],
                Polar=alignAndFocusArgs['Polar'],
                CalcInterfere=calc_interfere,
                SampleT=sample_t)

            ConvertToHistogram(
                InputWorkspace=sam_placzek,
                OutputWorkspace=sam_placzek)

        for wksp in [sam_placzek, sam_corrected]:
            ConvertUnits(
                InputWorkspace=wksp,
                OutputWorkspace=wksp,
                Target='MomentumTransfer',
                EMode='Elastic')

            Rebin(
                InputWorkspace=wksp,
                OutputWorkspace=wksp,
                Params=binning,
                PreserveEvents=True)

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

    # STEP 7:  S(Q) and F(Q), bank-by-bank

    sam_corrected_norm = sam_corrected + '_norm'
    to_absolute_scale(sam_corrected, sam_corrected_norm)
    save_banks(
        InputWorkspace=sam_corrected_norm,
        Filename=nexus_filename,
        Title="SQ_banks_normalized",
        OutputDir=OutputDir,
        GroupingWorkspace=grp_wksp,
        Binning=binning)

    if self_scattering_level_correction:
        sam_corrected_norm_scaled = sam_corrected_norm + '_scaled'
        _, bad_fitted_levels = \
            calculate_and_apply_fitted_levels(sam_corrected_norm,
                                              self_scattering_level_correction,
                                              sam_corrected_norm_scaled)
        if bad_fitted_levels:
            for bank, scale in bad_fitted_levels.items():
                final_message +=\
                    f'Bank {bank} potentially has a tilted baseline with ' \
                    f'the fitted scale = {scale} ' \
                    f'and was not scaled in {sam_corrected_norm_scaled}\n'
        save_banks(
            InputWorkspace=sam_corrected_norm_scaled,
            Filename=nexus_filename,
            Title="SQ_banks_normalized_scaled",
            OutputDir=OutputDir,
            GroupingWorkspace=grp_wksp,
            Binning=binning)
    else:
        sam_corrected_norm_scaled = sam_corrected_norm  # just an alias

    fq_banks = 'FQ_banks'
    to_f_of_q(sam_corrected_norm_scaled, fq_banks)
    save_banks(
        InputWorkspace=fq_banks,
        Filename=nexus_filename,
        Title="FQ_banks",
        OutputDir=OutputDir,
        GroupingWorkspace=grp_wksp,
        Binning=binning)

    # Print log information
    material = Material(sam_corrected)
    print("<b>^2:", material.bcoh_avg_sqrd)
    print("<b^2>:", material.btot_sqrd_avg)
    print("Laue term:", material.laue_monotonic_diffuse_scat)
    print("sample total xsection:", material.totalScatterXSection())
    print("vanadium total xsection:",
          Material(van_corrected).totalScatterXSection())

    # Output Bragg Diffraction
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
    xmax_rebin = max(xmax)
    tof_binning = "{xmin},-0.01,{xmax}".format(xmin=xmin_rebin, xmax=xmax_rebin)

    Rebin(
        InputWorkspace=sam_corrected,
        OutputWorkspace=sam_corrected,
        Params=tof_binning)

    SaveGSS(
        InputWorkspace=sam_corrected,
        Filename=os.path.join(os.path.abspath(OutputDir), title+".gsa"),
        SplitFiles=False,
        Append=False,
        MultiplyByBinWidth=True,
        Format="SLOG",
        ExtendedHeader=True)

    if final_message:
        log.warning(final_message)

    return mtd[sam_corrected]
