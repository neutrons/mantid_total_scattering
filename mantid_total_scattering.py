#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function)
import os
import six
import sys
import glob
import re
import json
import collections
from h5py import File
import mantid
from mantid import mtd
from mantid.simpleapi import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import m_n, micro
from scipy.constants import physical_constants
from scipy import interpolate, signal, ndimage, optimize

if six.PY3:
    unicode = str
    import configparser
else:
    import ConfigParser as configparser

from inelastic.placzek import CalculatePlaczekSelfScattering
from inelastic.placzek import GetIncidentSpectrumFromMonitor, FitIncidentSpectrum

#-----------------------------------------------------------------------------------------#
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

def save_file(ws, title, header=list()):
    with open(title, 'w') as f:
        for line in header:
            f.write('# %s \n' % line)
    SaveAscii(
        InputWorkspace=ws,
        Filename=title,
        Separator='Space',
        ColumnHeader=False,
        AppendToFile=True)


def save_banks_old(ws, title, binning=None):
    CloneWorkspace(InputWorkspace=ws, OutputWorkspace="tmp")
    # if mtd["tmp"].isDistribution():
    #    ConvertFromDistribution(mtd["tmp"])
    if binning:
        Rebin(InputWorkspace="tmp",
              OutputWorkspace="tmp",
              Params=binning,
              PreserveEvents=False)
    if mtd["tmp"].YUnit() == "Counts":
        try:
            print(
                "Unit:",
                mtd["tmp"].YUnit(),
                "Distribution:",
                mtd["tmp"].isDistribution())
            ConvertToDistribution("tmp")
        except BaseException:
            pass
    filename = os.path.join(os.getcwd(), title)
    print(filename)
    SaveAscii(InputWorkspace="tmp",
              Filename=filename,
              Separator='Space',
              ColumnHeader=False,
              AppendToFile=False,
              SpectrumList=range(mtd["tmp"].getNumberHistograms()))
    return


def save_banks(InputWorkspace, Filename, Title, OutputDir='./', Binning=None):
    CloneWorkspace(InputWorkspace=InputWorkspace, OutputWorkspace="tmp")
    # if mtd["tmp"].isDistribution():
    #    ConvertFromDistribution(mtd["tmp"])
    if Binning:
        Rebin(InputWorkspace="tmp",
              OutputWorkspace="tmp",
              Params=Binning,
              PreserveEvents=False)
    if mtd["tmp"].YUnit() == "Counts":
        try:
            print(
                "Unit:",
                mtd["tmp"].YUnit(),
                "Distribution:",
                mtd["tmp"].isDistribution())
            ConvertToDistribution("tmp")
        except BaseException:
            pass
    filename = os.path.join(OutputDir, Filename)
    print(filename)
    SaveNexusProcessed(
        InputWorkspace="tmp",
        Filename=Filename,
        Title=Title,
        Append=True,
        PreserveEvents=False,
        WorkspaceIndexList=range(
            mtd["tmp"].getNumberHistograms()))
    return


def save_banks_with_fit(
        title,
        fitrange_individual,
        InputWorkspace=None,
        **kwargs):
    # Header
    for i, fitrange in enumerate(fitrange_individual):
        print('fitrange:', fitrange[0], fitrange[1])

        Fit(Function='name=LinearBackground,A0=1.0,A1=0.0',
            WorkspaceIndex=i,
            # range cannot include area with NAN
            StartX=fitrange[0], EndX=fitrange[1],
            InputWorkspace=InputWorkspace, Output=InputWorkspace, OutputCompositeMembers=True)
        fitParams = mtd[InputWorkspace + '_Parameters']

        bank_title = title + '_' + InputWorkspace + '_bank_' + str(i) + '.dat'
        with open(bank_title, 'w') as f:
            if 'btot_sqrd_avg' in kwargs:
                f.write('#<b^2> : %f \n' % kwargs['btot_sqrd_avg'])
            if 'bcoh_avg_sqrd' in kwargs:
                f.write('#<b>^2 : %f \n' % kwargs['bcoh_avg_sqrd'])
            if 'self_scat' in kwargs:
                f.write('#self scattering : %f \n' % kwargs['self_scat'])
            f.write('#fitrange: %f %f \n' % (fitrange[0], fitrange[1]))
            f.write(
                '#for bank%d: %f + %f * Q\n' %
                (i +
                 1,
                 fitParams.cell(
                     'Value',
                     0),
                    fitParams.cell(
                     'Value',
                     1)))

    # Body
    for bank in range(mtd[InputWorkspace].getNumberHistograms()):
        x_data = mtd[InputWorkspace].readX(bank)[0:-1]
        y_data = mtd[InputWorkspace].readY(bank)
        bank_title = title + '_' + InputWorkspace + \
            '_bank_' + str(bank) + '.dat'
        print("####", bank_title)
        with open(bank_title, 'a') as f:
            for x, y in zip(x_data, y_data):
                f.write("%f %f \n" % (x, y))


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

#-------------------------------------------------------------------------
# Event Filters


def GenerateEventsFilterFromFiles(filenames, OutputWorkspace,
                                  InformationWorkspace, **kwargs):
    pass


def GenerateEventsFilterFromFiles(
        filenames,
        OutputWorkspace,
        InformationWorkspace,
        **kwargs):

    logName = kwargs.get('LogName', None)
    minValue = kwargs.get('MinimumLogValue', None)
    maxValue = kwargs.get('MaximumLogValue', None)
    logInterval = kwargs.get('LogValueInterval', None)
    unitOftime = kwargs.get('UnitOfTime', 'Nanoseconds')

    # TODO - handle multi-file filtering. Delete this line once implemented.
    assert len(
        filenames) == 1, 'ERROR: Multi-file filtering is not yet supported. (Stay tuned...)'

    for i, filename in enumerate(filenames):
        Load(Filename=filename, OutputWorkspace=filename)
        splitws, infows = GenerateEventsFilter(InputWorkspace=filename,
                                               UnitOfTime=unitOfTime,
                                               LogName=logName,
                                               MinimumLogValue=minValue,
                                               MaximumLogValue=maxValue,
                                               LogValueInterval=logInterval)
        if i == 0:
            GroupWorkspaces(splitws, OutputWorkspace=Outputworkspace)
            GroupWorkspaces(infows, OutputWorkspace=InformationWorkspace)
        else:
            mtd[OutputWorkspace].add(splitws)
            mtd[InformationWorkspace].add(infows)
    return

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
    if inelastic_dict is None:
        inelastic_dict = {"Type": None}
        return inelastic_dict

    corr_type = inelastic_dict["Type"]

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


#-------------------------------------------------------------------------
# MAIN - NOM_pdf


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Absolute normalization PDF generation")
    parser.add_argument('json', help='Input json file')
    options = parser.parse_args()

    print("loading config from '%s'" % options.json)
    with open(options.json, 'r') as handle:
        config = json.load(handle)
    title = config['Title']
    instr = config['Instrument']

    # Get sample info
    sample = config['Sample']
    sam_mass_density = sample.get('MassDensity', None)
    sam_packing_fraction = sample.get('PackingFraction', None)
    sam_geometry = sample.get('Geometry', None)
    sam_material = sample.get('Material', None)

    # Get normalization info
    van = config['Vanadium']
    van_mass_density = van.get('MassDensity', None)
    van_packing_fraction = van.get('PackingFraction', 1.0)
    van_geometry = van.get('Geometry', None)
    van_material = van.get('Material', 'V')

    # Get calibration, characterization, and other settings
    high_q_linear_fit_range = config['HighQLinearFitRange']
    binning = config['Merging']['QBinning']
    # workspace indices - zero indexed arrays
    wkspIndices = config['Merging']['SumBanks']
    # TODO how much of each bank gets merged has info here in the form of
    # {"ID", "Qmin", "QMax"}
    cache_dir = config.get("CacheDir", os.path.abspath('.'))
    output_dir = config.get("OutputDir", os.path.abspath('.'))

    # Create Nexus file basenames
    sample['Runs'] = procNumbers(sample['Runs'])
    sample['Background']['Runs'] = procNumbers(
        sample['Background'].get('Runs', None))

    sam_scans = ','.join(['%s_%d' % (instr, num) for num in sample['Runs']])
    container = ','.join(['%s_%d' % (instr, num)
                          for num in sample['Background']["Runs"]])
    container_bg = None
    if "Background" in sample['Background']:
        sample['Background']['Background']['Runs'] = procNumbers(
            sample['Background']['Background']['Runs'])
        container_bg = ','.join(['%s_%d' % (instr, num)
                                 for num in sample['Background']['Background']['Runs']])
        if len(container_bg) == 0:
            container_bg = None

    van['Runs'] = procNumbers(van['Runs'])
    van['Background']['Runs'] = procNumbers(van['Background']['Runs'])

    van_scans = ','.join(['%s_%d' % (instr, num) for num in van['Runs']])
    van_bg = ','.join(['%s_%d' % (instr, num)
                       for num in van['Background']["Runs"]])
    if len(van_bg) == 0:
        van_bg = None

    nexus_filename = title + '.nxs'
    try:
        os.remove(nexus_filename)
    except OSError:
        pass

    # Get absolute scale information from Nexus file
    print("#-----------------------------------#")
    print("# Sample")
    print("#-----------------------------------#")
    natoms, self_scat, sam_info = getAbsScaleInfoFromNexus(sample['Runs'],
                                                           PackingFraction=sam_packing_fraction,
                                                           SampleMassDensity=sam_mass_density,
                                                           Geometry=sam_geometry,
                                                           ChemicalFormula=sam_material)

    print("#-----------------------------------#")
    print("# Vanadium")
    print("#-----------------------------------#")
    van_geo_info = van_geometry.copy()
    nvan_atoms, tmp, van_info = getAbsScaleInfoFromNexus(van['Runs'],
                                                         PackingFraction=1.0,
                                                         SampleMassDensity=van_mass_density,
                                                         Geometry=van_geo_info,
                                                         ChemicalFormula="V")

    if natoms and nvan_atoms:
        print("Sample natoms:", natoms)
        print("Vanadium natoms:", nvan_atoms)
        print("Vanadium natoms / Sample natoms:", nvan_atoms / natoms)
        print()
    else:
        natoms = 1
        nvan_atoms
        print("No sample atoms or vanadium atoms => no normalization by atoms.")

    # Get sample corrections
    if sam_info:
        sam_mass_density = float(
            sample.get(
                'MassDensity',
                sam_info['mass_density']))
        sam_material = sample.get('Material', sam_info['formula'])
        sam_packing_fraction = sample.get(
            'PackingFraction', sam_info['packing_fraction'])
    sam_geometry = sample.get('Geometry', None)
    sam_abs_corr = sample.get("AbsorptionCorrection", None)
    sam_ms_corr = sample.get("MultipleScatteringCorrection", None)
    sam_inelastic_corr = SetInelasticCorrection(
        sample.get('InelasticCorrection', None))

    # Get vanadium corrections
    van_material = van.get('Material', 'V')
    van_mass_density = van.get('MassDensity', van_info['mass_density'])
    van_packing_fraction = van.get(
        'PackingFraction',
        van_info['packing_fraction'])
    van_geometry = van.get('Geometry', None)
    van_abs_corr = van.get("AbsorptionCorrection", {"Type": None})
    van_ms_corr = van.get("MultipleScatteringCorrection", {"Type": None})
    van_inelastic_corr = SetInelasticCorrection(
        van.get('InelasticCorrection', None))

    results = PDLoadCharacterizations(
        Filename=config['Merging']['Characterizations']['Filename'],
        OutputWorkspace='characterizations')
    alignAndFocusArgs = dict(PrimaryFlightPath=results[2],
                             SpectrumIDs=results[3],
                             L2=results[4],
                             Polar=results[5],
                             Azimuthal=results[6])

    alignAndFocusArgs['CalFilename'] = config['Calibration']['Filename']
    # alignAndFocusArgs['GroupFilename'] don't use
    # alignAndFocusArgs['Params'] use resampleX
    alignAndFocusArgs['ResampleX'] = -6000
    alignAndFocusArgs['Dspacing'] = True
    #alignAndFocusArgs['PreserveEvents'] = True
    alignAndFocusArgs['RemovePromptPulseWidth'] = 50
    alignAndFocusArgs['MaxChunkSize'] = 8
    # alignAndFocusArgs['CompressTolerance'] use defaults
    # alignAndFocusArgs['UnwrapRef'] POWGEN option
    # alignAndFocusArgs['LowResRef'] POWGEN option
    # alignAndFocusArgs['LowResSpectrumOffset'] POWGEN option
    # alignAndFocusArgs['CropWavelengthMin'] from characterizations file
    # alignAndFocusArgs['CropWavelengthMax'] from characterizations file
    alignAndFocusArgs['Characterizations'] = 'characterizations'
    alignAndFocusArgs['ReductionProperties'] = '__snspowderreduction'
    alignAndFocusArgs['CacheDir'] = os.path.abspath(cache_dir)

    # TODO take out the RecalculatePCharge in the future once tested

    #-----------------------------------------------------------------------------------------#
    # Load Sample
    AlignAndFocusPowderFromFiles(OutputWorkspace='sample',
                                 Filename=sam_scans,
                                 Absorption=None,
                                 **alignAndFocusArgs)
    sam_wksp = 'sample'
    NormaliseByCurrent(InputWorkspace=sam_wksp,
                       OutputWorkspace=sam_wksp,
                       RecalculatePCharge=True)

    ConvertUnits(InputWorkspace=sam_wksp,
                 OutputWorkspace=sam_wksp,
                 Target="MomentumTransfer",
                 EMode="Elastic")
    sample_title = "sample_and_container"
    print(os.path.join(output_dir, sample_title + ".dat"))
    save_banks(InputWorkspace=sam_wksp,
               Filename=nexus_filename,
               Title=sample_title,
               OutputDir=output_dir,
               Binning=binning)

    #-----------------------------------------------------------------------------------------#
    # Load Sample Container
    AlignAndFocusPowderFromFiles(OutputWorkspace='container',
                                 Filename=container,
                                 Absorption=None,
                                 **alignAndFocusArgs)
    container = 'container'
    NormaliseByCurrent(InputWorkspace=container,
                       OutputWorkspace=container,
                       RecalculatePCharge=True)
    ConvertUnits(InputWorkspace=container,
                 OutputWorkspace=container,
                 Target="MomentumTransfer",
                 EMode="Elastic")
    container_title = "container"
    save_banks(InputWorkspace=container,
               Filename=nexus_filename,
               Title=container_title,
               OutputDir=output_dir,
               Binning=binning)

    #-----------------------------------------------------------------------------------------#
    # Load Sample Container Background

    if container_bg is not None:
        AlignAndFocusPowderFromFiles(OutputWorkspace='container_background',
                                     Filename=container_bg,
                                     Absorption=None,
                                     **alignAndFocusArgs)
        container_bg = 'container_background'
        NormaliseByCurrent(InputWorkspace=container_bg,
                           OutputWorkspace=container_bg,
                           RecalculatePCharge=True)
        ConvertUnits(InputWorkspace=container_bg,
                     OutputWorkspace=container_bg,
                     Target="MomentumTransfer",
                     EMode="Elastic")
        container_bg_title = "container_background"
        save_banks(InputWorkspace=container_bg,
                   Filename=nexus_filename,
                   Title=container_bg_title,
                   OutputDir=output_dir,
                   Binning=binning)

    #-----------------------------------------------------------------------------------------#
    # Load Vanadium
    #Load(Filename=van_abs, OutputWorkspace='van_absorption')
    AlignAndFocusPowderFromFiles(OutputWorkspace='vanadium',
                                 Filename=van_scans,
                                 AbsorptionWorkspace=None,
                                 **alignAndFocusArgs)
    van_wksp = 'vanadium'
    if "Shape" not in van_geometry:
        van_geometry.update({'Shape': 'Cylinder'})
    van_geometry.update({'Center': [0., 0., 0., ]})
    NormaliseByCurrent(InputWorkspace=van_wksp,
                       OutputWorkspace=van_wksp,
                       RecalculatePCharge=True)
    SetSample(
        InputWorkspace=van_wksp,
        Geometry=van_geometry,
        Material={
            'ChemicalFormula': van_material,
            'SampleMassDensity': van_mass_density})
    ConvertUnits(InputWorkspace=van_wksp,
                 OutputWorkspace=van_wksp,
                 Target="MomentumTransfer",
                 EMode="Elastic")
    vanadium_title = "vanadium_and_background"

    save_banks(InputWorkspace=van_wksp,
               Filename=nexus_filename,
               Title=vanadium_title,
               OutputDir=output_dir,
               Binning=binning)

    #-----------------------------------------------------------------------------------------#
    # Load Vanadium Background
    if van_bg is not None:
        AlignAndFocusPowderFromFiles(OutputWorkspace='vanadium_background',
                                     Filename=van_bg,
                                     AbsorptionWorkspace=None,
                                     **alignAndFocusArgs)

    van_bg = 'vanadium_background'
    NormaliseByCurrent(InputWorkspace=van_bg,
                       OutputWorkspace=van_bg,
                       RecalculatePCharge=True)
    ConvertUnits(InputWorkspace=van_bg,
                 OutputWorkspace=van_bg,
                 Target="MomentumTransfer",
                 EMode="Elastic")
    vanadium_bg_title = "vanadium_background"
    save_banks(InputWorkspace=van_bg,
               Filename=nexus_filename,
               Title=vanadium_bg_title,
               OutputDir=output_dir,
               Binning=binning)

    #-----------------------------------------------------------------------------------------#
    # Load Instrument Characterizations
    PDDetermineCharacterizations(InputWorkspace=sam_wksp,
                                 Characterizations='characterizations',
                                 ReductionProperties='__snspowderreduction')
    propMan = PropertyManagerDataService.retrieve('__snspowderreduction')
    qmax = 2. * np.pi / propMan['d_min'].value
    qmin = 2. * np.pi / propMan['d_max'].value
    for a, b in zip(qmin, qmax):
        print('Qrange:', a, b)
    mask_info = generateCropingTable(qmin, qmax)

    #-----------------------------------------------------------------------------------------#
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
               OutputDir=output_dir,
               Binning=binning)
    save_banks(InputWorkspace=van_wksp,
               Filename=nexus_filename,
               Title=vanadium_title,
               OutputDir=output_dir,
               Binning=binning)
    save_banks(InputWorkspace=sam_wksp,
               Filename=nexus_filename,
               Title=sample_title,
               OutputDir=output_dir,
               Binning=binning)

    #-----------------------------------------------------------------------------------------#
    # STEP 2.0: Prepare vanadium as normalization calibrant

    # Multiple-Scattering and Absorption (Steps 2-4) for Vanadium

    van_corrected = 'van_corrected'
    ConvertUnits(InputWorkspace=van_wksp,
                 OutputWorkspace=van_corrected,
                 Target="Wavelength",
                 EMode="Elastic")

    if "Type" in van_abs_corr:
        if van_abs_corr['Type'] == 'Carpenter' or van_ms_corr['Type'] == 'Carpenter':
            MultipleScatteringCylinderAbsorption(
                InputWorkspace=van_corrected,
                OutputWorkspace=van_corrected,
                CylinderSampleRadius=van['Geometry']['Radius'])
        elif van_abs_corr['Type'] == 'Mayers' or van_ms_corr['Type'] == 'Mayers':
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
               OutputDir=output_dir,
               Binning=binning)
    save_banks(InputWorkspace=van_corrected,
               Filename=nexus_filename,
               Title=vanadium_title + "_with_peaks",
               OutputDir=output_dir,
               Binning=binning)

    # TODO subtract self-scattering of vanadium (According to Eq. 7 of Howe,
    # McGreevey, and Howells, JPCM, 1989)

    # Smooth Vanadium (strip peaks plus smooth)

    ConvertUnits(InputWorkspace=van_corrected,
                 OutputWorkspace=van_corrected,
                 Target='dSpacing',
                 EMode='Elastic')
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
               OutputDir=output_dir,
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
               OutputDir=output_dir,
               Binning=binning)

    # Inelastic correction
    if van_inelastic_corr['Type'] == "Placzek":
        for van_scan in van['Runs']:
            van_incident_wksp = 'van_incident_wksp'
            lambda_binning_fit = van['InelasticCorrection']['LambdaBinningForFit']
            lambda_binning_calc = van['InelasticCorrection']['LambdaBinningForCalc']
            GetIncidentSpectrumFromMonitor(
                '%s_%s' %
                (instr, str(van_scan)), OutputWorkspace=van_incident_wksp)

            fit_type = van['InelasticCorrection']['FitSpectrumWith']
            FitIncidentSpectrum(InputWorkspace=van_incident_wksp,
                                OutputWorkspace=van_incident_wksp,
                                FitSpectrumWith=fit_type,
                                BinningForFit=lambda_binning_fit,
                                BinningForCalc=lambda_binning_calc,
                                plot_diagnostics=False)

            van_placzek = 'van_placzek'

            SetSample(InputWorkspace=van_incident_wksp,
                      Material={'ChemicalFormula': van_material,
                                'SampleMassDensity': van_mass_density})
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
                   OutputDir=output_dir,
                   Binning=binning)

        # Rebin in Wavelength
        for wksp in [van_placzek, van_corrected]:
            ConvertUnits(InputWorkspace=wksp,
                         OutputWorkspace=wksp,
                         Target='Wavelength',
                         EMode='Elastic')
            Rebin(InputWorkspace=wksp, OutputWorkspace=wksp,
                  Params=lambda_binning_calc, PreserveEvents=False)

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
                   OutputDir=output_dir,
                   Binning=binning)

    ConvertUnits(InputWorkspace=van_corrected,
                 OutputWorkspace=van_corrected,
                 Target='MomentumTransfer',
                 EMode='Elastic')

    SetUncertainties(InputWorkspace=van_corrected,
                     OutputWorkspace=van_corrected,
                     SetError='zero')

    #-----------------------------------------------------------------------------------------#
    # STEP 2.1: Normalize by Vanadium

    for name in [sam_wksp, van_corrected]:
        ConvertUnits(
            InputWorkspace=name,
            OutputWorkspace=name,
            Target='MomentumTransfer',
            EMode='Elastic',
            ConvertFromPointData=False)
        Rebin(InputWorkspace=name, OutputWorkspace=name,
              Params=binning, PreserveEvents=False)
        if not mtd[name].isDistribution():
            ConvertToDistribution(name)

    Divide(
        LHSWorkspace=sam_wksp,
        RHSWorkspace=van_corrected,
        OutputWorkspace=sam_wksp)
    Divide(
        LHSWorkspace=sam_raw,
        RHSWorkspace=van_corrected,
        OutputWorkspace=sam_raw)

    sample_title += "_normalized"
    save_banks(InputWorkspace=sam_wksp,
               Filename=nexus_filename,
               Title=sample_title,
               OutputDir=output_dir,
               Binning=binning)
    save_banks(InputWorkspace=sam_raw,
               Filename=nexus_filename,
               Title="sample_normalized",
               OutputDir=output_dir,
               Binning=binning)

    for name in [container, van_corrected]:
        ConvertUnits(
            InputWorkspace=name,
            OutputWorkspace=name,
            Target='MomentumTransfer',
            EMode='Elastic',
            ConvertFromPointData=False)
        Rebin(InputWorkspace=name, OutputWorkspace=name,
              Params=binning, PreserveEvents=False)
        if not mtd[name].isDistribution():
            ConvertToDistribution(name)
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
    Divide(
        LHSWorkspace=container_raw,
        RHSWorkspace=van_corrected,
        OutputWorkspace=container_raw)
    if van_bg is not None:
        Divide(
            LHSWorkspace=van_bg,
            RHSWorkspace=van_corrected,
            OutputWorkspace=van_bg)
    if container_bg is not None:
        Divide(
            LHSWorkspace=container_bg,
            RHSWorkspace=van_corrected,
            OutputWorkspace=container_bg)

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
               OutputDir=output_dir,
               Binning=binning)
    save_banks(InputWorkspace=container_raw,
               Filename=nexus_filename,
               Title="container_normalized",
               OutputDir=output_dir,
               Binning=binning)

    if container_bg is not None:
        container_bg_title += "_normalised"
        save_banks(InputWorkspace=container_bg,
                   Filename=nexus_filename,
                   Title=container_bg_title,
                   OutputDir=output_dir,
                   Binning=binning)

    if van_bg is not None:
        vanadium_bg_title += "_normalized"
        save_banks(InputWorkspace=van_bg,
                   Filename=nexus_filename,
                   Title=vanadium_bg_title,
                   OutputDir=output_dir,
                   Binning=binning)

    #-----------------------------------------------------------------------------------------#
    # STEP 3 & 4: Subtract multiple scattering and apply absorption correction

    ConvertUnits(
        InputWorkspace=sam_wksp,
        OutputWorkspace=sam_wksp,
        Target="Wavelength",
        EMode="Elastic")

    sam_corrected = 'sam_corrected'
    if sam_abs_corr:
        if sam_abs_corr['Type'] == 'Carpenter' or sam_ms_corr['Type'] == 'Carpenter':
            MultipleScatteringCylinderAbsorption(
                InputWorkspace=sam_wksp,
                OutputWorkspace=sam_corrected,
                CylinderSampleRadius=sample['Geometry']['Radius'])
        elif sam_abs_corr['Type'] == 'Mayers' or sam_ms_corr['Type'] == 'Mayers':
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
                   OutputDir=output_dir,
                   Binning=binning)
    else:
        CloneWorkspace(InputWorkspace=sam_wksp, OutputWorkspace=sam_corrected)

    #-----------------------------------------------------------------------------------------#
    # STEP 5: Divide by number of atoms in sample

    mtd[sam_corrected] = (nvan_atoms / natoms) * mtd[sam_corrected]
    ConvertUnits(InputWorkspace=sam_corrected, OutputWorkspace=sam_corrected,
                 Target='MomentumTransfer', EMode='Elastic')
    sample_title += "_norm_by_atoms"
    save_banks(InputWorkspace=sam_corrected,
               Filename=nexus_filename,
               Title=sample_title,
               OutputDir=output_dir,
               Binning=binning)

    #-----------------------------------------------------------------------------------------#
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
               OutputDir=output_dir,
               Binning=binning)

    #-----------------------------------------------------------------------------------------#
    # STEP 7: Inelastic correction
    ConvertUnits(InputWorkspace=sam_corrected, OutputWorkspace=sam_corrected,
                 Target='Wavelength', EMode='Elastic')
    if sam_inelastic_corr['Type'] == "Placzek":
        if sam_material is None:
            raise Exception(
                "ERROR: For Placzek correction, must specifiy a sample material.")
        for sam_scan in sample['Runs']:
            sam_incident_wksp = 'sam_incident_wksp'
            lambda_binning_fit = sample['InelasticCorrection']['LambdaBinningForFit']
            lambda_binning_calc = sample['InelasticCorrection']['LambdaBinningForCalc']
            GetIncidentSpectrumFromMonitor(
                '%s_%s' %
                (instr, str(sam_scan)), OutputWorkspace=sam_incident_wksp)

            fit_type = sample['InelasticCorrection']['FitSpectrumWith']
            FitIncidentSpectrum(InputWorkspace=sam_incident_wksp,
                                OutputWorkspace=sam_incident_wksp,
                                FitSpectrumWith=fit_type,
                                BinningForFit=lambda_binning_fit,
                                BinningForCalc=lambda_binning_calc)

            sam_placzek = 'sam_placzek'
            SetSample(InputWorkspace=sam_incident_wksp,
                      Material={'ChemicalFormula': sam_material,
                                'SampleMassDensity': sam_mass_density})
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
                   OutputDir=output_dir,
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
                   OutputDir=output_dir,
                   Binning=binning)

    #-----------------------------------------------------------------------------------------#
    # STEP 7: Output spectrum

    # TODO Since we already went from Event -> 2D workspace, can't use this anymore
    # if alignAndFocusArgs['PreserveEvents']:
    #    CompressEvents(InputWorkspace=sam_corrected, OutputWorkspace=sam_corrected)
    #    CompressEvents(InputWorkspace=van_corrected, OutputWorkspace=van_corrected)

    #-----------------------------------------------------------------------------------------#

    # F(Q) bank-by-bank Section
    CloneWorkspace(InputWorkspace=sam_corrected, OutputWorkspace='FQ_banks_ws')
    FQ_banks = 'FQ_banks'

    # S(Q) bank-by-bank Section
    material = mtd[sam_corrected].sample().getMaterial()
    if material.name() is None or len(material.name().strip()) == 0:
        raise RuntimeError('Sample material was not set')
    bcoh_avg_sqrd = material.cohScatterLength() * material.cohScatterLength()
    btot_sqrd_avg = material.totalScatterLengthSqrd()
    laue_monotonic_diffuse_scat = btot_sqrd_avg / bcoh_avg_sqrd
    CloneWorkspace(InputWorkspace=sam_corrected, OutputWorkspace='SQ_banks_ws')
    SQ_banks = (1. / bcoh_avg_sqrd) * \
        mtd['SQ_banks_ws'] - laue_monotonic_diffuse_scat + 1.

    save_banks(InputWorkspace="FQ_banks_ws",
               Filename=nexus_filename,
               Title="FQ_banks",
               OutputDir=output_dir,
               Binning=binning)
    save_banks(InputWorkspace="SQ_banks_ws",
               Filename=nexus_filename,
               Title="SQ_banks",
               OutputDir=output_dir,
               Binning=binning)

    #-----------------------------------------------------------------------------------------#
    # STOP HERE FOR NOW
    print("<b>^2:", bcoh_avg_sqrd)
    print("<b^2>:", btot_sqrd_avg)
    print("Laue term:", laue_monotonic_diffuse_scat)
    print(mtd[sam_corrected].sample().getMaterial().totalScatterXSection())
    print(mtd[van_corrected].sample().getMaterial().totalScatterXSection())

    '''
    #-----------------------------------------------------------------------------------------#
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

    save_banks('SQ_banks',     title=os.path.join(output_dir,title+"_SQ_banks.dat"),     binning=binning)
    save_banks('FQ_banks',     title=os.path.join(output_dir,title+"_FQ_banks.dat"),     binning=binning)
    save_banks('FQ_banks_raw', title=os.path.join(output_dir,title+"_FQ_banks_raw.dat"), binning=binning)

    #-----------------------------------------------------------------------------------------#
    # Event workspace -> Histograms
    Rebin(InputWorkspace=sam_corrected, OutputWorkspace=sam_corrected, Params=binning, PreserveEvents=True)
    Rebin(InputWorkspace=van_corrected, OutputWorkspace=van_corrected, Params=binning, PreserveEvents=True)
    Rebin(InputWorkspace='container',   OutputWorkspace='container',   Params=binning, PreserveEvents=True)
    Rebin(InputWorkspace='sample',      OutputWorkspace='sample',      Params=binning, PreserveEvents=True)
    if van_bg is not None:
        Rebin(InputWorkspace=van_bg,        OutputWorkspace='background',      Params=binning, PreserveEvents=True)

    #-----------------------------------------------------------------------------------------#
    # Apply Qmin Qmax limits

    #MaskBinsFromTable(InputWorkspace=sam_corrected, OutputWorkspace='sam_single',       MaskingInformation=mask_info)
    #MaskBinsFromTable(InputWorkspace=van_corrected, OutputWorkspace='van_single',       MaskingInformation=mask_info)
    #MaskBinsFromTable(InputWorkspace='container',   OutputWorkspace='container_single', MaskingInformation=mask_info)
    #MaskBinsFromTable(InputWorkspace='sample',      OutputWorkspace='sample_raw_single',MaskingInformation=mask_info)

    #-----------------------------------------------------------------------------------------#
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

    #-----------------------------------------------------------------------------------------#
    # Merged S(Q) and F(Q)

    save_banks(InputWorkspace="FQ_banks_ws",
               Filename=nexus_filename,
               Title="FQ_banks",
               OutputDir=output_dir,
               Binning=binning)
    save_banks(InputWorkspace="SQ_banks_ws",
               Filename=nexus_filename,
               Title="SQ_banks",
               OutputDir=output_dir,
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
                    'for merged banks %s: %f + %f * Q' % (','.join([ str(i) for i in wkspIndices]), \
                                                       fitParams.cell('Value', 0), fitParams.cell('Value', 1)) ]
'''
