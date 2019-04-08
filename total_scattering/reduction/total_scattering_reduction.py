#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function)

import os
import itertools
import numpy as np
from scipy.constants import Avogadro

from mantid import mtd
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
    SetSample, \
    SetUncertainties, \
    StripVanadiumPeaks

from total_scattering.file_handling.load import load
from total_scattering.file_handling.save import save_banks
from total_scattering.inelastic.placzek import \
    CalculatePlaczekSelfScattering, \
    FitIncidentSpectrum, \
    GetIncidentSpectrumFromMonitor

# Utilities


def myMatchingBins(leftWorkspace, rightWorkspace):
    leftXData = mtd[leftWorkspace].dataX(0)
    rightXData = mtd[rightWorkspace].dataX(0)

    if len(leftXData) != len(rightXData):
        return False

    if abs(sum(leftXData) - sum(rightXData)) > 1.e-7:
        print(
            "Sums do not match: LHS = ",
            sum(leftXData),
            "RHS =",
            sum(rightXData))
        return False

    leftDeltaX = leftXData[0] - leftXData[1]
    rightDeltaX = rightXData[0] - rightXData[1]

    if abs(leftDeltaX -
           rightDeltaX) >= 1e-4 or abs(rightXData[0] -
                                       leftXData[0]) >= 1e-4:
        return False

    return True


def generateCropingTable(qmin, qmax):
    mask_info = CreateEmptyTableWorkspace()
    mask_info.addColumn("str", "SpectraList")
    mask_info.addColumn("double", "XMin")
    mask_info.addColumn("double", "XMax")
    for (i, value) in enumerate(qmin):
        mask_info.addRow([str(i), 0.0, value])
    for (i, value) in enumerate(qmax):
        mask_info.addRow([str(i), value, 100.0])

    return mask_info


def getQmaxFromData(Workspace=None, WorkspaceIndex=0):
    if Workspace is None:
        return None
    return max(mtd[Workspace].readX(WorkspaceIndex))

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
                                "Interference": False,
                                "FitSpectrumWith": "GaussConvCubicSpline",
                                "LambdaBinning": "0.16,0.04,2.8"}
            inelastic_settings = default_settings.copy()
            inelastic_settings.update(inelastic_dict)

        else:
            raise Exception("Unknown Inelastic Correction Type")

    return inelastic_settings


def TotalScatteringReduction(config=None):
    facility = config['Facility']
    title = config['Title']
    instr = config['Instrument']

    # Get sample info
    sample = config['Sample']
    sam_mass_density = sample.get('MassDensity', None)
    sam_packing_fraction = sample.get('PackingFraction', None)
    sam_geometry = sample.get('Geometry', None)
    sam_material = sample.get('Material', None)

    # Get normalization info
    van = config['Normalization']
    van_mass_density = van.get('MassDensity', None)
    van_packing_fraction = van.get('PackingFraction', 1.0)
    van_geometry = van.get('Geometry', None)
    van_material = van.get('Material', 'V')

    # Get calibration, characterization, and other settings
    merging = config['Merging']
    binning = merging['QBinning']

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

    # Get vanadium corrections
    van_mass_density = van.get('MassDensity', van_mass_density)
    van_packing_fraction = van.get(
        'PackingFraction',
        van_packing_fraction)
    van_abs_corr = van.get("AbsorptionCorrection", {"Type": None})
    van_ms_corr = van.get("MultipleScatteringCorrection", {"Type": None})
    van_inelastic_corr = SetInelasticCorrection(
        van.get('InelasticCorrection', None))

    alignAndFocusArgs = dict()
    alignAndFocusArgs['CalFilename'] = config['Calibration']['Filename']
    # alignAndFocusArgs['GroupFilename'] don't use
    # alignAndFocusArgs['Params'] = "0.,0.02,40."
    alignAndFocusArgs['ResampleX'] = -6000
    alignAndFocusArgs['Dspacing'] = False
    alignAndFocusArgs['PreserveEvents'] = False
    alignAndFocusArgs['MaxChunkSize'] = 8
    alignAndFocusArgs['CacheDir'] = os.path.abspath(cache_dir)

    # Get any additional AlignAndFocusArgs from JSON input
    if "AlignAndFocusArgs" in config:
        otherArgs = config["AlignAndFocusArgs"]
        alignAndFocusArgs.update(otherArgs)

    print(alignAndFocusArgs)
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
        **alignAndFocusArgs)
    sample_title = "sample_and_container"
    print(os.path.join(OutputDir, sample_title + ".dat"))
    print("HERE:", mtd[sam_wksp].getNumberHistograms())
    print(grp_wksp)
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
    container = load('container', container_scans, **alignAndFocusArgs)
    save_banks(InputWorkspace=container,
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
        save_banks(InputWorkspace=container_bg,
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
        **alignAndFocusArgs)
    vanadium_title = "vanadium_and_background"

    save_banks(InputWorkspace=van_wksp,
               Filename=nexus_filename,
               Title=vanadium_title,
               OutputDir=OutputDir,
               GroupingWorkspace=grp_wksp,
               Binning=binning)

    van_molecular_mass = mtd[van_wksp].sample(
    ).getMaterial().relativeMolecularMass()
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
        save_banks(InputWorkspace=van_bg,
                   Filename=nexus_filename,
                   Title=vanadium_bg_title,
                   OutputDir=OutputDir,
                   GroupingWorkspace=grp_wksp,
                   Binning=binning)

    # Load Instrument Characterizations
    if "Characterizations" in merging:
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
        # mask_info = generateCropingTable(qmin, qmax)

    # STEP 1: Subtract Backgrounds

    sam_raw = 'sam_raw'
    CloneWorkspace(
        InputWorkspace=sam_wksp,
        OutputWorkspace=sam_raw)  # for later

    container_raw = 'container_raw'
    CloneWorkspace(InputWorkspace=container,
                   OutputWorkspace=container_raw)  # for later

    if van_bg is not None:
        Minus(
            LHSWorkspace=van_wksp,
            RHSWorkspace=van_bg,
            OutputWorkspace=van_wksp)
    Minus(
        LHSWorkspace=sam_wksp,
        RHSWorkspace=container,
        OutputWorkspace=sam_wksp)
    if container_bg is not None:
        Minus(
            LHSWorkspace=container,
            RHSWorkspace=container_bg,
            OutputWorkspace=container)

    for wksp in [container, van_wksp, sam_wksp]:
        ConvertUnits(InputWorkspace=wksp,
                     OutputWorkspace=wksp,
                     Target="MomentumTransfer",
                     EMode="Elastic")
    container_title = "container_minus_back"
    vanadium_title = "vanadium_minus_back"
    sample_title = "sample_minus_back"
    save_banks(InputWorkspace=container,
               Filename=nexus_filename,
               Title=container_title,
               OutputDir=OutputDir,
               GroupingWorkspace=grp_wksp,
               Binning=binning)
    save_banks(InputWorkspace=van_wksp,
               Filename=nexus_filename,
               Title=vanadium_title,
               OutputDir=OutputDir,
               GroupingWorkspace=grp_wksp,
               Binning=binning)
    save_banks(InputWorkspace=sam_wksp,
               Filename=nexus_filename,
               Title=sample_title,
               OutputDir=OutputDir,
               GroupingWorkspace=grp_wksp,
               Binning=binning)

    # STEP 2.0: Prepare vanadium as normalization calibrant

    # Multiple-Scattering and Absorption (Steps 2-4) for Vanadium

    van_corrected = 'van_corrected'
    ConvertUnits(InputWorkspace=van_wksp,
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
                MayersSampleCorrection(InputWorkspace=van_corrected,
                                       OutputWorkspace=van_corrected,
                                       MultipleScattering=True)
            else:
                MayersSampleCorrection(InputWorkspace=van_corrected,
                                       OutputWorkspace=van_corrected,
                                       MultipleScattering=False)
        else:
            print("NO VANADIUM absorption or multiple scattering!")
    else:
        CloneWorkspace(
            InputWorkspace=van_corrected,
            OutputWorkspace=van_corrected)

    ConvertUnits(InputWorkspace=van_corrected,
                 OutputWorkspace=van_corrected,
                 Target='MomentumTransfer',
                 EMode='Elastic')
    vanadium_title += "_ms_abs_corrected"
    save_banks(InputWorkspace=van_corrected,
               Filename=nexus_filename,
               Title=vanadium_title,
               OutputDir=OutputDir,
               GroupingWorkspace=grp_wksp,
               Binning=binning)
    save_banks(InputWorkspace=van_corrected,
               Filename=nexus_filename,
               Title=vanadium_title + "_with_peaks",
               OutputDir=OutputDir,
               GroupingWorkspace=grp_wksp,
               Binning=binning)

    # TODO subtract self-scattering of vanadium (According to Eq. 7 of Howe,
    # McGreevey, and Howells, JPCM, 1989)

    # Smooth Vanadium (strip peaks plus smooth)

    ConvertUnits(InputWorkspace=van_corrected,
                 OutputWorkspace=van_corrected,
                 Target='dSpacing',
                 EMode='Elastic')
    # After StripVanadiumPeaks, the workspace goes from EventWorkspace ->
    # Workspace2D
    StripVanadiumPeaks(InputWorkspace=van_corrected,
                       OutputWorkspace=van_corrected,
                       BackgroundType='Quadratic')
    ConvertUnits(InputWorkspace=van_corrected,
                 OutputWorkspace=van_corrected,
                 Target='MomentumTransfer',
                 EMode='Elastic')
    vanadium_title += '_peaks_stripped'
    save_banks(InputWorkspace=van_corrected,
               Filename=nexus_filename,
               Title=vanadium_title,
               OutputDir=OutputDir,
               GroupingWorkspace=grp_wksp,
               Binning=binning)

    ConvertUnits(InputWorkspace=van_corrected,
                 OutputWorkspace=van_corrected,
                 Target='TOF',
                 EMode='Elastic')
    FFTSmooth(InputWorkspace=van_corrected,
              OutputWorkspace=van_corrected,
              Filter="Butterworth",
              Params='20,2',
              IgnoreXBins=True,
              AllSpectra=True)
    ConvertUnits(InputWorkspace=van_corrected,
                 OutputWorkspace=van_corrected,
                 Target='MomentumTransfer',
                 EMode='Elastic')
    vanadium_title += '_smoothed'
    save_banks(InputWorkspace=van_corrected,
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
        FitIncidentSpectrum(InputWorkspace=van_incident_wksp,
                            OutputWorkspace=van_incident_wksp,
                            FitSpectrumWith=fit_type,
                            BinningForFit=lambda_binning_fit,
                            BinningForCalc=lambda_binning_calc,
                            PlotDiagnostics=False)

        van_placzek = 'van_placzek'

        SetSample(InputWorkspace=van_incident_wksp,
                  Material={'ChemicalFormula': str(van_material),
                            'SampleMassDensity': str(van_mass_density)})
        CalculatePlaczekSelfScattering(IncidentWorkspace=van_incident_wksp,
                                       ParentWorkspace=van_corrected,
                                       OutputWorkspace=van_placzek,
                                       L1=19.5,
                                       L2=alignAndFocusArgs['L2'],
                                       Polar=alignAndFocusArgs['Polar'])
        ConvertToHistogram(InputWorkspace=van_placzek,
                           OutputWorkspace=van_placzek)

        # Save before rebin in Q
        for wksp in [van_placzek, van_corrected]:
            ConvertUnits(InputWorkspace=wksp,
                         OutputWorkspace=wksp,
                         Target='MomentumTransfer',
                         EMode='Elastic')
            Rebin(InputWorkspace=wksp, OutputWorkspace=wksp,
                  Params=binning, PreserveEvents=True)

        save_banks(InputWorkspace=van_placzek,
                   Filename=nexus_filename,
                   Title="vanadium_placzek",
                   OutputDir=OutputDir,
                   GroupingWorkspace=grp_wksp,
                   Binning=binning)

        # Rebin in Wavelength
        for wksp in [van_placzek, van_corrected]:
            ConvertUnits(InputWorkspace=wksp,
                         OutputWorkspace=wksp,
                         Target='Wavelength',
                         EMode='Elastic')
            Rebin(InputWorkspace=wksp, OutputWorkspace=wksp,
                  Params=lambda_binning_calc, PreserveEvents=True)

        # Save after rebin in Q
        for wksp in [van_placzek, van_corrected]:
            ConvertUnits(InputWorkspace=wksp,
                         OutputWorkspace=wksp,
                         Target='MomentumTransfer',
                         EMode='Elastic')

        # Subtract correction in Wavelength
        for wksp in [van_placzek, van_corrected]:
            ConvertUnits(InputWorkspace=wksp,
                         OutputWorkspace=wksp,
                         Target='Wavelength',
                         EMode='Elastic')
            if not mtd[wksp].isDistribution():
                ConvertToDistribution(wksp)

        Minus(LHSWorkspace=van_corrected,
              RHSWorkspace=van_placzek,
              OutputWorkspace=van_corrected)

        # Save after subtraction
        for wksp in [van_placzek, van_corrected]:
            ConvertUnits(InputWorkspace=wksp,
                         OutputWorkspace=wksp,
                         Target='MomentumTransfer',
                         EMode='Elastic')
        vanadium_title += '_placzek_corrected'
        save_banks(InputWorkspace=van_corrected,
                   Filename=nexus_filename,
                   Title=vanadium_title,
                   OutputDir=OutputDir,
                   GroupingWorkspace=grp_wksp,
                   Binning=binning)

    ConvertUnits(InputWorkspace=van_corrected,
                 OutputWorkspace=van_corrected,
                 Target='MomentumTransfer',
                 EMode='Elastic')

    SetUncertainties(InputWorkspace=van_corrected,
                     OutputWorkspace=van_corrected,
                     SetError='zero')

    # STEP 2.1: Normalize by Vanadium

    for name in [sam_wksp, van_corrected]:
        ConvertUnits(
            InputWorkspace=name,
            OutputWorkspace=name,
            Target='MomentumTransfer',
            EMode='Elastic',
            ConvertFromPointData=False)
        Rebin(InputWorkspace=name, OutputWorkspace=name,
              Params=binning, PreserveEvents=True)
        # if not mtd[name].isDistribution():
        #    ConvertToDistribution(name)

    Divide(
        LHSWorkspace=sam_wksp,
        RHSWorkspace=van_corrected,
        OutputWorkspace=sam_wksp)
    # Divide(
    #    LHSWorkspace=sam_raw,
    #    RHSWorkspace=van_corrected,
    #    OutputWorkspace=sam_raw)

    sample_title += "_normalized"
    save_banks(InputWorkspace=sam_wksp,
               Filename=nexus_filename,
               Title=sample_title,
               OutputDir=OutputDir,
               GroupingWorkspace=grp_wksp,
               Binning=binning)
    save_banks(InputWorkspace=sam_raw,
               Filename=nexus_filename,
               Title="sample_normalized",
               OutputDir=OutputDir,
               GroupingWorkspace=grp_wksp,
               Binning=binning)

    for name in [container, van_corrected]:
        ConvertUnits(
            InputWorkspace=name,
            OutputWorkspace=name,
            Target='MomentumTransfer',
            EMode='Elastic',
            ConvertFromPointData=False)
        Rebin(InputWorkspace=name, OutputWorkspace=name,
              Params=binning, PreserveEvents=True)
        # if not mtd[name].isDistribution():
        #    ConvertToDistribution(name)
    print()
    print("## Container ##")
    print("YUnit:", mtd[container].YUnit(), "|", mtd[van_corrected].YUnit())
    print(
        "blocksize:",
        mtd[container].blocksize(),
        mtd[van_corrected].blocksize())
    print("dist:", mtd[container].isDistribution(),
          mtd[van_corrected].isDistribution())
    print("Do bins match?:", myMatchingBins(container, van_corrected))
    print(
        "Distributions?",
        mtd[container].isDistribution(),
        mtd[van_corrected].isDistribution())
    print()

    Divide(
        LHSWorkspace=container,
        RHSWorkspace=van_corrected,
        OutputWorkspace=container)
    # Divide(
    #    LHSWorkspace=container_raw,
    #    RHSWorkspace=van_corrected,
    #    OutputWorkspace=container_raw)
    # if van_bg is not None:
    #    Divide(
    #        LHSWorkspace=van_bg,
    #        RHSWorkspace=van_corrected,
    #        OutputWorkspace=van_bg)
    # if container_bg is not None:
    #    Divide(
    #       LHSWorkspace=container_bg,
    #        RHSWorkspace=van_corrected,
    #        OutputWorkspace=container_bg)

    print()
    print("## Container After Divide##")
    print("YUnit:", mtd[container].YUnit(), "|", mtd[van_corrected].YUnit())
    print(
        "blocksize:",
        mtd[container].blocksize(),
        mtd[van_corrected].blocksize())
    print("dist:", mtd[container].isDistribution(),
          mtd[van_corrected].isDistribution())
    print("Do bins match?:", myMatchingBins(container, van_corrected))
    print(
        "Distributions?",
        mtd[container].isDistribution(),
        mtd[van_corrected].isDistribution())
    print()

    container_title += '_normalized'
    save_banks(InputWorkspace=container,
               Filename=nexus_filename,
               Title=container_title,
               OutputDir=OutputDir,
               GroupingWorkspace=grp_wksp,
               Binning=binning)
    save_banks(InputWorkspace=container_raw,
               Filename=nexus_filename,
               Title="container_normalized",
               OutputDir=OutputDir,
               GroupingWorkspace=grp_wksp,
               Binning=binning)

    # if container_bg is not None:
    #    container_bg_title += "_normalised"
    #    save_banks(InputWorkspace=container_bg,
    #               Filename=nexus_filename,
    #               Title=container_bg_title,
    #               OutputDir=OutputDir,
    #               GroupingWorkspace=grp_wksp,
    #               Binning=binning)

    # if van_bg is not None:
    #    vanadium_bg_title += "_normalized"
    #    save_banks(InputWorkspace=van_bg,
    #               Filename=nexus_filename,
    #               Title=vanadium_bg_title,
    #               OutputDir=OutputDir,
    #               GroupingWorkspace=grp_wksp,
    #               Binning=binning)

    # STEP 3 & 4: Subtract multiple scattering and apply absorption correction

    ConvertUnits(
        InputWorkspace=sam_wksp,
        OutputWorkspace=sam_wksp,
        Target="Wavelength",
        EMode="Elastic")

    sam_corrected = 'sam_corrected'
    if sam_abs_corr:
        if sam_abs_corr['Type'] == 'Carpenter' \
                or sam_ms_corr['Type'] == 'Carpenter':
            CarpenterSampleCorrection(
                InputWorkspace=sam_wksp,
                OutputWorkspace=sam_corrected,
                CylinderSampleRadius=sample['Geometry']['Radius'])
        elif sam_abs_corr['Type'] == 'Mayers' \
                or sam_ms_corr['Type'] == 'Mayers':
            if sam_ms_corr['Type'] == 'Mayers':
                MayersSampleCorrection(InputWorkspace=sam_wksp,
                                       OutputWorkspace=sam_corrected,
                                       MultipleScattering=True)
            else:
                MayersSampleCorrection(InputWorkspace=sam_wksp,
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
        save_banks(InputWorkspace=sam_corrected,
                   Filename=nexus_filename,
                   Title=sample_title,
                   OutputDir=OutputDir,
                   GroupingWorkspace=grp_wksp,
                   Binning=binning)
    else:
        CloneWorkspace(InputWorkspace=sam_wksp, OutputWorkspace=sam_corrected)

    # STEP 5: Divide by number of atoms in sample

    mtd[sam_corrected] = (nvan_atoms / natoms) * mtd[sam_corrected]
    ConvertUnits(InputWorkspace=sam_corrected, OutputWorkspace=sam_corrected,
                 Target='MomentumTransfer', EMode='Elastic')
    sample_title += "_norm_by_atoms"
    save_banks(InputWorkspace=sam_corrected,
               Filename=nexus_filename,
               Title=sample_title,
               OutputDir=OutputDir,
               GroupingWorkspace=grp_wksp,
               Binning=binning)

    # STEP 6: Divide by total scattering length squared = total scattering
    # cross-section over 4 * pi
    sigma_v = mtd[van_corrected].sample().getMaterial().totalScatterXSection()
    prefactor = (sigma_v / (4. * np.pi))
    print("Total scattering cross-section of Vanadium:",
          sigma_v, " sigma_v / 4*pi:", prefactor)
    mtd[sam_corrected] = prefactor * mtd[sam_corrected]
    sample_title += '_multiply_by_vanSelfScat'
    save_banks(InputWorkspace=sam_corrected,
               Filename=nexus_filename,
               Title=sample_title,
               OutputDir=OutputDir,
               GroupingWorkspace=grp_wksp,
               Binning=binning)

    # STEP 7: Inelastic correction
    ConvertUnits(InputWorkspace=sam_corrected, OutputWorkspace=sam_corrected,
                 Target='Wavelength', EMode='Elastic')
    if sam_inelastic_corr['Type'] == "Placzek":
        if sam_material is None:
            error = "For Placzek correction, must specifiy a sample material."
            raise Exception(error)
        for sam_scan in sample['Runs']:
            sam_incident_wksp = 'sam_incident_wksp'
            sam_inelastic_opts = sample['InelasticCorrection']
            lambda_binning_fit = sam_inelastic_opts['LambdaBinningForFit']
            lambda_binning_calc = sam_inelastic_opts['LambdaBinningForCalc']
            GetIncidentSpectrumFromMonitor(
                Filename=facility_file_format % (instr, sam_scan),
                OutputWorkspace=sam_incident_wksp)

            fit_type = sample['InelasticCorrection']['FitSpectrumWith']
            FitIncidentSpectrum(InputWorkspace=sam_incident_wksp,
                                OutputWorkspace=sam_incident_wksp,
                                FitSpectrumWith=fit_type,
                                BinningForFit=lambda_binning_fit,
                                BinningForCalc=lambda_binning_calc)

            sam_placzek = 'sam_placzek'
            SetSample(InputWorkspace=sam_incident_wksp,
                      Material={'ChemicalFormula': str(sam_material),
                                'SampleMassDensity': str(sam_mass_density)})
            CalculatePlaczekSelfScattering(IncidentWorkspace=sam_incident_wksp,
                                           ParentWorkspace=sam_corrected,
                                           OutputWorkspace=sam_placzek,
                                           L1=19.5,
                                           L2=alignAndFocusArgs['L2'],
                                           Polar=alignAndFocusArgs['Polar'])
            ConvertToHistogram(InputWorkspace=sam_placzek,
                               OutputWorkspace=sam_placzek)

        # Save before rebin in Q
        for wksp in [sam_placzek, sam_corrected]:
            ConvertUnits(InputWorkspace=wksp,
                         OutputWorkspace=wksp,
                         Target='MomentumTransfer',
                         EMode='Elastic')
            Rebin(InputWorkspace=wksp, OutputWorkspace=wksp,
                  Params=binning, PreserveEvents=True)

        save_banks(InputWorkspace=sam_placzek,
                   Filename=nexus_filename,
                   Title="sample_placzek",
                   OutputDir=OutputDir,
                   GroupingWorkspace=grp_wksp,
                   Binning=binning)

        # Save after rebin in Q
        for wksp in [sam_placzek, sam_corrected]:
            ConvertUnits(InputWorkspace=wksp,
                         OutputWorkspace=wksp,
                         Target='MomentumTransfer',
                         EMode='Elastic')

        Minus(LHSWorkspace=sam_corrected,
              RHSWorkspace=sam_placzek,
              OutputWorkspace=sam_corrected)

        # Save after subtraction
        for wksp in [sam_placzek, sam_corrected]:
            ConvertUnits(InputWorkspace=wksp,
                         OutputWorkspace=wksp,
                         Target='MomentumTransfer',
                         EMode='Elastic')
        sample_title += '_placzek_corrected'
        save_banks(InputWorkspace=sam_corrected,
                   Filename=nexus_filename,
                   Title=sample_title,
                   OutputDir=OutputDir,
                   GroupingWorkspace=grp_wksp,
                   Binning=binning)

    # STEP 7: Output spectrum

    # TODO Since we already went from Event -> 2D workspace, can't use this
    # anymore
    print('sam:', mtd[sam_corrected].id())
    print('van:', mtd[van_corrected].id())
    if alignAndFocusArgs['PreserveEvents']:
        CompressEvents(
            InputWorkspace=sam_corrected,
            OutputWorkspace=sam_corrected)

    # F(Q) bank-by-bank Section
    CloneWorkspace(InputWorkspace=sam_corrected, OutputWorkspace='FQ_banks_ws')
    # TODO: Add the following when implemented - FQ_banks = 'FQ_banks'

    # S(Q) bank-by-bank Section
    material = mtd[sam_corrected].sample().getMaterial()
    if material.name() is None or len(material.name().strip()) == 0:
        raise RuntimeError('Sample material was not set')
    bcoh_avg_sqrd = material.cohScatterLength() * material.cohScatterLength()
    btot_sqrd_avg = material.totalScatterLengthSqrd()
    laue_monotonic_diffuse_scat = btot_sqrd_avg / bcoh_avg_sqrd
    CloneWorkspace(InputWorkspace=sam_corrected, OutputWorkspace='SQ_banks_ws')

    # TODO: Add the following when implemented
    '''
    SQ_banks = (1. / bcoh_avg_sqrd) * \
        mtd['SQ_banks_ws'] - laue_monotonic_diffuse_scat + 1.
    '''

    save_banks(InputWorkspace="FQ_banks_ws",
               Filename=nexus_filename,
               Title="FQ_banks",
               OutputDir=OutputDir,
               GroupingWorkspace=grp_wksp,
               Binning=binning)
    save_banks(InputWorkspace="SQ_banks_ws",
               Filename=nexus_filename,
               Title="SQ_banks",
               OutputDir=OutputDir,
               GroupingWorkspace=grp_wksp,
               Binning=binning)

    # STOP HERE FOR NOW
    print("<b>^2:", bcoh_avg_sqrd)
    print("<b^2>:", btot_sqrd_avg)
    print("Laue term:", laue_monotonic_diffuse_scat)
    print(
        "sample total xsection:",
        mtd[sam_corrected].sample().getMaterial().totalScatterXSection())
    print(
        "vanadium total xsection:",
        mtd[van_corrected].sample().getMaterial().totalScatterXSection())

    # Output Bragg Diffraction
    ConvertUnits(InputWorkspace=sam_corrected,
                 OutputWorkspace=sam_corrected,
                 Target="TOF",
                 EMode="Elastic")

    ConvertToHistogram(InputWorkspace=sam_corrected,
                       OutputWorkspace=sam_corrected)

    Rebin(InputWorkspace=sam_corrected,
          OutputWorkspace=sam_corrected,
          Params="350.0,-0.0001,26233.0")
    xmin = "450.0,1100.0,2050.0,3400.0,4500.0"
    xmax = "11000.0,11000.0,11000.0,11000.0,11000.0"
    CropWorkspaceRagged(InputWorkspace=sam_corrected,
                        OutputWorkspace=sam_corrected,
                        Xmin=xmin,
                        Xmax=xmax)
    return mtd[sam_corrected]
    # ResampleX(InputWorkspace=sam_corrected,
    #          OutputWorkspace=sam_corrected,
    #          NumberBins=3000,
    #          LogBinning=True)

    # SaveGSS(InputWorkspace=sam_corrected,
    #        Filename=os.path.join(OutputDir,title+".gsa"),
    #        SplitFiles=False,
    #        Append=False,
    #        MultiplyByBinWidth=True,
    #        Format="SLOG",
    #        ExtendedHeader=True)
    # process the run
    '''
    SNSPowderReduction(
        Filename=sam_scans,
        MaxChunkSize=alignAndFocusArgs['MaxChunkSize'],
        PreserveEvents=True,
        PushDataPositive='ResetToZero',
        CalibrationFile=alignAndFocusArgs['CalFilename'],
        CharacterizationRunsFile=merging['Characterizations']['Filename'],
        BackgroundNumber=sample["Background"]["Runs"],
        VanadiumNumber=van["Runs"],
        VanadiumBackgroundNumber=van["Background"]["Runs"],
        RemovePromptPulseWidth=alignAndFocusArgs['RemovePromptPulseWidth'],
        ResampleX=alignAndFocusArgs['ResampleX'],
        BinInDspace=True,
        FilterBadPulses=25.,
        SaveAs="gsas fullprof topas",
        OutputFilePrefix=title,
        OutputDirectory=OutputDir,
        StripVanadiumPeaks=True,
        VanadiumRadius=van_geometry['Radius'],
        NormalizeByCurrent=True,
        FinalDataUnits="dSpacing")

    # Ouput bank-by-bank with linear fits for high-Q

    # fit the last 80% of the bank being used
    for i, q in zip(range(mtd[sam_corrected].getNumberHistograms()), qmax):
        qmax_data = getQmaxFromData(sam_corrected, i)
        qmax[i] = q if q <= qmax_data else qmax_data

    fitrange_individual = [(high_q_linear_fit_range*q, q) for q in qmax]

    for q in qmax:
        print('Linear Fit Qrange:', high_q_linear_fit_range*q, q)


    kwargs = { 'btot_sqrd_avg' : btot_sqrd_avg,
               'bcoh_avg_sqrd' : bcoh_avg_sqrd,
               'self_scat' : self_scat }

    save_banks_with_fit( title, fitrange_individual, InputWorkspace='SQ_banks', **kwargs)
    save_banks_with_fit( title, fitrange_individual, InputWorkspace='FQ_banks', **kwargs)
    save_banks_with_fit( title, fitrange_individual, InputWorkspace='FQ_banks_raw', **kwargs)

    save_banks('SQ_banks', title=os.path.join(OutputDir,title+"_SQ_banks.dat"),
                binning=binning)
    save_banks('FQ_banks', title=os.path.join(OutputDir,title+"_FQ_banks.dat"),
                binning=binning)
    save_banks('FQ_banks_raw', title=os.path.join(OutputDir,title+"_FQ_banks_raw.dat"),
                binning=binning)

    # Event workspace -> Histograms
    Rebin(InputWorkspace=sam_corrected, OutputWorkspace=sam_corrected,
          Params=binning, PreserveEvents=True)
    Rebin(InputWorkspace=van_corrected, OutputWorkspace=van_corrected,
          Params=binning, PreserveEvents=True)
    Rebin(InputWorkspace='container',   OutputWorkspace='container',
          Params=binning, PreserveEvents=True)
    Rebin(InputWorkspace='sample', OutputWorkspace='sample',
          Params=binning, PreserveEvents=True)
    if van_bg is not None:
        Rebin(InputWorkspace=van_bg, OutputWorkspace='background',
              Params=binning, PreserveEvents=True)

    # Apply Qmin Qmax limits

    #MaskBinsFromTable(InputWorkspace=sam_corrected, OutputWorkspace='sam_single',
                       MaskingInformation=mask_info)
    #MaskBinsFromTable(InputWorkspace=van_corrected, OutputWorkspace='van_single',
                       MaskingInformation=mask_info)
    #MaskBinsFromTable(InputWorkspace='container', OutputWorkspace='container_single',
                       MaskingInformation=mask_info)
    #MaskBinsFromTable(InputWorkspace='sample', OutputWorkspace='sample_raw_single',
                       MaskingInformation=mask_info)

    # Get sinlge, merged spectrum from banks

    CloneWorkspace(InputWorkspace=sam_corrected, OutputWorkspace='sam_single')
    CloneWorkspace(InputWorkspace=van_corrected, OutputWorkspace='van_single')
    CloneWorkspace(InputWorkspace='container', OutputWorkspace='container_single')
    CloneWorkspace(InputWorkspace='sample', OutputWorkspace='sample_raw_single')
    CloneWorkspace(InputWorkspace='background', OutputWorkspace='background_single')

    SumSpectra(InputWorkspace='sam_single', OutputWorkspace='sam_single',
               ListOfWorkspaceIndices=wkspIndices)
    SumSpectra(InputWorkspace='van_single', OutputWorkspace='van_single',
               ListOfWorkspaceIndices=wkspIndices)

    # Diagnostic workspaces
    SumSpectra(InputWorkspace='container_single', OutputWorkspace='container_single',
               ListOfWorkspaceIndices=wkspIndices)
    SumSpectra(InputWorkspace='sample_raw_single', OutputWorkspace='sample_raw_single',
               ListOfWorkspaceIndices=wkspIndices)
    SumSpectra(InputWorkspace='background_single', OutputWorkspace='background_single',
               ListOfWorkspaceIndices=wkspIndices)

    # Merged S(Q) and F(Q)

    save_banks(InputWorkspace="FQ_banks_ws",
               Filename=nexus_filename,
               Title="FQ_banks",
               OutputDir=OutputDir,
               Binning=binning)
    save_banks(InputWorkspace="SQ_banks_ws",
               Filename=nexus_filename,
               Title="SQ_banks",
               OutputDir=OutputDir,
               Binning=binning)


    # do the division correctly and subtract off the material specific term
    CloneWorkspace(InputWorkspace='sam_single', OutputWorkspace='SQ_ws')
    SQ = (1./bcoh_avg_sqrd)*mtd['SQ_ws'] - (term_to_subtract-1.)  # +1 to get back to S(Q)

    CloneWorkspace(InputWorkspace='sam_single', OutputWorkspace='FQ_ws')
    FQ_raw = mtd['FQ_ws']
    FQ = FQ_raw - self_scat

    qmax = 48.0
    Fit(Function='name=LinearBackground,A0=1.0,A1=0.0',
        StartX=high_q_linear_fit_range*qmax, EndX=qmax, # range cannot include area with NAN
        InputWorkspace='SQ', Output='SQ', OutputCompositeMembers=True)
    fitParams = mtd['SQ_Parameters']

    qmax = getQmaxFromData('FQ', WorkspaceIndex=0)
    Fit(Function='name=LinearBackground,A0=1.0,A1=0.0',
        StartX=high_q_linear_fit_range*qmax, EndX=qmax, # range cannot include area with NAN
        InputWorkspace='FQ', Output='FQ', OutputCompositeMembers=True)
    fitParams = mtd['FQ_Parameters']

    qmax = 48.0
    Fit(Function='name=LinearBackground,A0=1.0,A1=0.0',
        StartX=high_q_linear_fit_range*qmax, EndX=qmax, # range cannot include area with NAN
        InputWorkspace='FQ_raw', Output='FQ_raw', OutputCompositeMembers=True)
    fitParams = mtd['FQ_raw_Parameters']

    # Save dat file
    header_lines = ['<b^2> : %f ' % btot_sqrd_avg, \
                    '<b>^2 : %f ' % bcoh_avg_sqrd, \
                    'self scattering: %f ' % self_scat, \
                    'fitrange: %f %f '  % (high_q_linear_fit_range*qmax,qmax), \
                    'for merged banks %s: %f + %f * Q' \
                    % (','.join([ str(i) for i in wkspIndices]), \
                    fitParams.cell('Value', 0), fitParams.cell('Value', 1)) ]
'''
