import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, ndimage, interpolate, optimize

from mantid import mtd
from mantid.simpleapi import \
    CalculateEfficiencyCorrection, \
    ConvertToPointData, \
    ConvertUnits, \
    CreateWorkspace, \
    Divide, \
    LoadNexusMonitors, \
    Multiply, \
    Rebin, \
    ResampleX, \
    SplineSmoothing

# ----------------------------------------------------------------------------------------- #
# Functions for fitting the incident spectrum


def getFitRange(x, y, x_lo, x_hi):
    if x_lo is None:
        x_lo = min(x)
    if x_hi is None:
        x_hi = max(x)

    x_fit = x[(x >= x_lo) & (x <= x_hi)]
    y_fit = y[(x >= x_lo) & (x <= x_hi)]
    return x_fit, y_fit


def fitCubicSpline(x_fit, y_fit, x, s=1e15):
    tck = interpolate.splrep(x_fit, y_fit, s=s)
    fit = interpolate.splev(x, tck, der=0)
    fit_prime = interpolate.splev(x, tck, der=1)
    return fit, fit_prime


def fitCubicSplineViaMantidSplineSmoothing(InputWorkspace, Params, **kwargs):
    Rebin(
        InputWorkspace=InputWorkspace,
        OutputWorkspace='fit',
        Params=Params,
        PreserveEvents=True)
    SplineSmoothing(
        InputWorkspace='fit',
        OutputWorkspace='fit',
        OutputWorkspaceDeriv='fit_prime',
        DerivOrder=1,
        **kwargs)
    return mtd['fit'].readY(0), mtd['fit_prime_1'].readY(0)


def fitHowellsFunction(x_fit, y_fit, x):
    # Fit with analytical function from HowellsEtAl
    def calc_HowellsFunction(
            lambdas,
            phi_max,
            phi_epi,
            lam_t,
            lam_1,
            lam_2,
            a):
        term1 = phi_max * ((lam_t**4.) / lambdas**5.) * \
            np.exp(-(lam_t / lambdas)**2.)
        term2 = (phi_epi / (lambdas**(1. + 2. * a))) * \
            (1. / (1 + np.exp((lambdas - lam_1) / lam_2)))
        return term1 + term2

    def calc_HowellsFunction1stDerivative(
            lambdas, phi_max, phi_epi, lam_t, lam_1, lam_2, a):
        term1 = (((2 * lam_t**2) / lambdas**2) - 5.) * (1. / lambdas) * \
            phi_max * ((lam_t**4.) / lambdas**5.) * np.exp(-(lam_t / lambdas)**2.)
        term2 = ((1 + 2 * a) / lambdas) \
            * (1. / lambdas) * (phi_epi / (lambdas ** (1. + 2. * a))) \
            * (1. / (1 + np.exp((lambdas - lam_1) / lam_2)))
        return term1 + term2

    params = [1., 1., 1., 0., 1., 1.]
    params, convergence = optimize.curve_fit(
        calc_HowellsFunction, x_fit, y_fit, params)
    fit = calc_HowellsFunction(x, *params)
    fit_prime = calc_HowellsFunction1stDerivative(x, *params)
    return fit, fit_prime


def fitCubicSplineWithGaussConv(x_fit, y_fit, x, sigma=3):
    # Fit with Cubic Spline using a Gaussian Convolution to get weights
    def moving_average(y, sigma=sigma):
        b = signal.gaussian(39, sigma)
        average = ndimage.filters.convolve1d(y, b / b.sum())
        var = ndimage.filters.convolve1d(np.power(y - average, 2), b / b.sum())
        return average, var

    avg, var = moving_average(y_fit)
    spline_fit = interpolate.UnivariateSpline(x_fit, y_fit, w=1. / np.sqrt(var))
    spline_fit_prime = spline_fit.derivative()
    fit = spline_fit(x)
    fit_prime = spline_fit_prime(x)
    return fit, fit_prime


# ----------------------------------------------------------------------------------------- #
# Get incident spectrum from Monitor

def plotIncidentSpectrum(x, y, x_fit, fit, fit_prime, title=None):
    plt.plot(x, y, 'bo', x_fit, fit, '--')
    plt.legend(['Incident Spectrum', 'Fit f(x)'], loc='best')
    if title is not None:
        plt.title(title)
    plt.show()

    plt.plot(x_fit, fit_prime / fit, 'x--', label="Fit f'(x)/f(x)")
    plt.xlabel('Wavelength')
    plt.legend()
    if title is not None:
        plt.title(title)
    axes = plt.gca()
    axes.set_ylim([-12, 6])
    plt.show()
    return


def GetIncidentSpectrumFromMonitor(
        Filename,
        OutputWorkspace="IncidentWorkspace",
        IncidentIndex=0,
        TransmissionIndex=1,
        Binning=".1,6000,2.9",
        BinType="ResampleX"):

    # -------------------------------------------------
    # Joerg's read_bm.pro code

    # Loop workspaces to get each incident spectrum
    monitor = 'monitor'
    LoadNexusMonitors(Filename=Filename, OutputWorkspace=monitor)
    ConvertUnits(InputWorkspace=monitor, OutputWorkspace=monitor,
                 Target='Wavelength', EMode='Elastic')
    lambdaMin, lambdaBinning, lambdaMax = [float(x) for x in Binning.split(',')]
    for x in [lambdaMin, lambdaBinning, lambdaMax]:
        print(x, type(x))
    if BinType == 'ResampleX':
        ResampleX(InputWorkspace=monitor,
                  OutputWorkspace=monitor,
                  XMin=[lambdaMin],  # TODO change ResampleX
                  XMax=[lambdaMax],
                  NumberBins=abs(int(lambdaBinning)),
                  LogBinning=(int(lambdaBinning) < 0),
                  PreserveEvents=True)
    elif BinType == 'Rebin':
        Rebin(InputWorkspace=monitor,
              OutputWorkspace=monitor,
              Params=[lambdaMin, lambdaBinning, lambdaMax],
              PreserveEvents=True)
    ConvertToPointData(InputWorkspace=monitor, OutputWorkspace=monitor)

    lam = mtd[monitor].readX(IncidentIndex)    # wavelength in A
    bm = mtd[monitor].readY(IncidentIndex)     # neutron counts / microsecond
    p = 0.000794807                       # Pressure (empirically adjusted to match eff.)
    thickness = .1                        # 1 mm = .1 cm
    abs_xs_3He = 5333.0                   # barns for lambda == 1.798 A
    p_to_rho = 2.43e-5                    # pressure to rho (atoms/angstroms^3)
    # p is set to give efficiency of 1.03 10^-5 at 1.8 A
    e0 = abs_xs_3He * lam / 1.798 * p_to_rho * p * thickness
    print('Efficiency:', 1. - np.exp(-e0))
    bmeff = bm / (1. - np.exp(-e0))      # neutron counts / microsecond
    print(bmeff)
    # bmeff = bmeff / constants.micro      # neutron counts / second

    CreateWorkspace(DataX=lam, DataY=bmeff,
                    OutputWorkspace=OutputWorkspace, UnitX='Wavelength')
    mtd[OutputWorkspace].setYUnit('Counts')
    return mtd[OutputWorkspace]


def FitIncidentSpectrum(InputWorkspace, OutputWorkspace,
                        FitSpectrumWith='GaussConvCubicSpline',
                        BinningForFit="0.15,0.05,3.2",
                        BinningForCalc=None,
                        PlotDiagnostics=False):

    incident_ws = mtd[InputWorkspace]

    # Fit Incident Spectrum
    # Get axis for actual calc (either provided in BinningForCalc or extracted
    # from incident wksp)
    incident_index = 0
    if BinningForCalc is None:
        x = incident_ws.readX(incident_index)
    else:
        try:
            params = [float(x) for x in BinningForCalc.split(',')]
        except AttributeError:
            params = [float(x) for x in BinningForCalc]
        xlo, binsize, xhi = params
        x = np.arange(xlo, xhi, binsize)

    Rebin(
        incident_ws,
        OutputWorkspace='fit',
        Params=BinningForFit,
        PreserveEvents=True)
    x_fit = np.array(mtd['fit'].readX(incident_index))
    y_fit = np.array(mtd['fit'].readY(incident_index))

    if FitSpectrumWith == 'CubicSpline':
        fit, fit_prime = fitCubicSpline(x_fit, y_fit, x, s=1e7)
        if PlotDiagnostics:
            plotIncidentSpectrum(
                x_fit,
                y_fit,
                x,
                fit,
                fit_prime,
                title='Simple Cubic Spline: Default')

    elif FitSpectrumWith == 'CubicSplineViaMantid':
        fit, fit_prime = fitCubicSplineViaMantidSplineSmoothing(
            InputWorkspace, Params=BinningForFit, MaxNumberOfBreaks=8)
        if PlotDiagnostics:
            plotIncidentSpectrum(
                x_fit,
                y_fit,
                x,
                fit,
                fit_prime,
                title='Cubic Spline via Mantid SplineSmoothing')

    elif FitSpectrumWith == 'HowellsFunction':
        fit, fit_prime = fitHowellsFunction(x_fit, y_fit, x)
        if PlotDiagnostics:
            plotIncidentSpectrum(
                x_fit,
                y_fit,
                x,
                fit,
                fit_prime,
                title='HowellsFunction')

    elif FitSpectrumWith == 'GaussConvCubicSpline':
        fit, fit_prime = fitCubicSplineWithGaussConv(x_fit, y_fit, x, sigma=2)
        if PlotDiagnostics:
            plotIncidentSpectrum(
                x_fit,
                y_fit,
                x,
                fit,
                fit_prime,
                title='Cubic Spline w/ Gaussian Kernel Convolution ')

    else:
        raise Exception("Unknown method for fitting incident spectrum")
        return

    CreateWorkspace(
        DataX=x,
        DataY=np.append(
            fit,
            fit_prime),
        OutputWorkspace=OutputWorkspace,
        UnitX='Wavelength',
        NSpec=2,
        Distribution=False)
    return mtd[OutputWorkspace]


if '__main__' == __name__:
    run_file = False
    nomad_test = True
    if run_file:
        # ------------------------------------------------------------------------------------- #
        # Get input parameters
        configfile = sys.argv[1]
        with open(configfile) as handle:
            config = json.loads(handle.read())

        sample = config['Sample']
        opts = sample['InelasticCorrection']

        # ------------------------------------------------------------------------------------- #
        # Get incident spectrum test for NOMAD
        runs = sample["Runs"].split(',')
        runs = ["%s_%s" % (config["Instrument"], run) for run in runs]

        binning = "0.05,0.00022,3.5"

        if nomad_test:
            runs[0] = "NOM_33943"
            binning = "0.0212406,0.00022,3.39828"  # matches read_bm.pro for lambda[100:15999]
        print("Processing Scan: ", runs[0])

        monitor = 'monitor'
        incident_ws = 'incident_ws'

        fig, (ax_bm, ax_bmeff) = plt.subplots(2, subplot_kw={'projection': 'mantid'}, sharex=True)

        # Beam Monitor
        LoadNexusMonitors(Filename=runs[0], OutputWorkspace=monitor)
        ConvertUnits(InputWorkspace=monitor, OutputWorkspace=monitor,
                     Target='Wavelength', EMode='Elastic')

        Rebin(InputWorkspace=monitor,
              OutputWorkspace=monitor,
              Params=binning,
              PreserveEvents=False)

        ax_bm.plot(mtd[monitor], '-', wkspIndex=0, label='Monitor', distribution=True)
        ax_bm.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

        # Use sample info
        CalculateEfficiencyCorrection(InputWorkspace=monitor,
                                      ChemicalFormula="(He3)",
                                      DensityType="Number Density",
                                      Density=1.93138101e-08,
                                      Thickness=.1,
                                      OutputWorkspace=incident_ws)

        Divide(LHSWorkspace=monitor, RHSWorkspace=incident_ws, OutputWorkspace=incident_ws)

        ax_bmeff.plot(mtd[incident_ws], '-', wkspIndex=0, label='Incident Spectrum (density)',
                      distribution=True)

        # Use measured efficiency
        CalculateEfficiencyCorrection(InputWorkspace=monitor,
                                      ChemicalFormula="(He3)",
                                      Efficiency=1.03e-5,
                                      OutputWorkspace=incident_ws)

        Divide(LHSWorkspace=monitor, RHSWorkspace=incident_ws, OutputWorkspace=incident_ws)

        ax_bmeff.plot(mtd[incident_ws], 'o', wkspIndex=0, label='Incident Spectrum (efficiency)',
                      distribution=True)

        # Use alpha
        CalculateEfficiencyCorrection(InputWorkspace=monitor,
                                      Alpha=-5.72861786781e-06,
                                      OutputWorkspace=incident_ws)

        Divide(LHSWorkspace=monitor, RHSWorkspace=incident_ws, OutputWorkspace=incident_ws)

        ax_bmeff.plot(mtd[incident_ws], '--', wkspIndex=0, label='Incident Spectrum (alpha)',
                      distribution=True)

        # Plot all
        ax_bm.legend()
        ax_bmeff.legend()
        plt.show()
        exit()
    # ----------------------------------------------------------------------------------------- #
    # Howells incident spectrum for testing
    howells_test = True
    if howells_test:
        # Howells function
        def delta(lam, lam_1, lam_2, lam_3=None, lam_4=None, alpha=None):
            retVal = 1. / (1. + np.exp((lam - lam_1) / lam_2))

            if lam_3 and lam_4 and alpha:
                term_2 = 1. + (alpha / (1. + np.exp((lam_3 - lam) / lam_4)))
                retVal *= term_2

            return retVal

        def phi_m(lam, **kwargs):
            phi_max = kwargs['phi_max']
            lam_t = kwargs['lam_t']
            return phi_max * (lam_t**4. / lam**5.) * np.exp(-(lam_t / lam)**2.)

        def phi_e(lam, **kwargs):
            phi_epi = kwargs['phi_epi']
            alpha = kwargs['alpha']
            lam_1 = kwargs['lam_1']
            lam_2 = kwargs['lam_2']

            a = None
            lam_3 = None
            lam_4 = None
            if 'a' in kwargs:
                a = kwargs['a']
            if 'lam_3' in kwargs:
                lam_3 = kwargs['lam_3']
            if 'lam_4' in kwargs:
                lam_4 = kwargs['lam_4']

            if lam_3 and lam_4:
                delta_term = delta(lam, lam_1, lam_2, lam_3, lam_4, a)
            else:
                delta_term = delta(lam, lam_1, lam_2)
            return phi_epi * delta_term / (lam**(1 + 2 * alpha))

        def calc_HowellsFunction(lam, **kwargs):
            return phi_m(lam, **kwargs) + phi_e(lam, **kwargs)

        incident_spectrums = dict()
        incident_spectrums['Ambient 300K polyethylene'] = {
            'phi_max': 6324,
            'phi_epi': 786,
            'lam_t': 1.58,
            'alpha': 0.099,
            'lam_1': 0.67143,
            'lam_2': 0.06075
        }

        incident_spectrums['Ambient 300K poisoned (Gd "meat" in polyethylene slabs)'] = {
            'phi_max': 1200,
            'phi_epi': 786,
            'lam_t': 1.58,
            'alpha': 0.099,
            'lam_1': 0.67143,
            'lam_2': 0.06075
        }

        incident_spectrums['Cold 77K polyethylene (Howells)'] = {
            'phi_max': 3838,
            'phi_epi': 1029,
            'lam_t': 2.97,
            'alpha': 0.089,
            'lam_1': 1.3287,
            'lam_2': 0.14735
        }

        incident_spectrums['Cold 77K polyethylene (MildnerEtAl)'] = {
            'phi_max': 3838,
            'phi_epi': 1029,
            'lam_t': 2.97,
            'alpha': 0.089,
            'lam_1': 1.3287,
            'lam_2': 0.14735,
            'a': 0.2882,
            'lam_3': 0.7253,
            'lam_4': 0.0486
        }

        fig, (ax_bm, ax_eff) = plt.subplots(2, subplot_kw={'projection': 'mantid'}, sharex=True)

        # Using Mantid
        lam_lo = 0.2
        lam_hi = 4.0
        lam_delta = 0.01
        binning = "%s,%s,%s" % (lam_lo, lam_delta, lam_hi)
        for moderator, spectrum_params in incident_spectrums.items():
            color = next(ax_bm._get_lines.prop_cycler)['color']

            incident_ws = 'howells_%s' % moderator
            CreateWorkspace(OutputWorkspace=incident_ws, NSpec=1, DataX=[0], DataY=[0],
                            UnitX='Wavelength', VerticalAxisUnit='Text',
                            VerticalAxisValues='IncidentSpectrum')
            Rebin(InputWorkspace=incident_ws, OutputWorkspace=incident_ws, Params=binning)
            ConvertToPointData(InputWorkspace=incident_ws, OutputWorkspace=incident_ws)

            wavelengths = mtd[incident_ws].readX(0)
            incident_spectrum = calc_HowellsFunction(wavelengths, **spectrum_params)
            mtd[incident_ws].setY(0, incident_spectrum)
            ax_bm.plot(mtd[incident_ws], '-', color=color, wkspIndex=0, label=moderator)

            eff_ws = 'efficiency'
            CalculateEfficiencyCorrection(InputWorkspace=incident_ws,
                                          Alpha=-0.693,
                                          OutputWorkspace=eff_ws)
            ConvertToPointData(InputWorkspace=eff_ws, OutputWorkspace=eff_ws)

            ax_eff.plot(mtd[eff_ws], '-', color=color, wkspIndex=0, label=moderator + ' efficiency')

            sample_ws = 'sample_ws'
            Multiply(LHSWorkspace=incident_ws, RHSWorkspace=eff_ws, OutputWorkspace=sample_ws)

            ax_bm.plot(mtd[sample_ws], 'o',
                       color=color,
                       wkspIndex=0,
                       label=moderator + ' measurement')

        ax_bm.legend()
        ax_eff.legend()
        plt.show()

    exit()
    # ----------------------------------------------------------------------------------------- #
    # Fit incident spectrum
    incident_fit = 'incident_fit'
    fit_type = sample['InelasticCorrection']['FitSpectrumWith']
    FitIncidentSpectrum(InputWorkspace=incident_ws,
                        OutputWorkspace=incident_fit,
                        FitSpectrumWith=fit_type,
                        BinningForFit=opts['LambdaBinningForFit'],
                        BinningForCalc=opts['LambdaBinningForCalc'],
                        PlotDiagnostics=opts['PlotFittingDiagnostics'])
