import sys
import json
import collections
import numpy as np
import scipy
from mantid.simpleapi import *

from inelastic.incident_spectrum import get_incident_spectrum_from_monitor, fit_incident_spectrum


# -------------------------------------------------------------------------
# Placzek - 1st order inelastic correction

def plot_placzek(x, y, fit, fit_prime, title=None):
    plt.plot(x, y, 'bo', x, fit, '--')
    plt.legend(['Incident Spectrum', 'Fit f(x)'], loc='best')
    if title is not None:
        plt.title(title)
    plt.show()

    plt.plot(x, x * fit_prime / fit, 'x--', label="Fit x*f'(x)/f(x)")
    plt.xlabel('Wavelength')
    plt.legend()
    if title is not None:
        plt.title(title)
    plt.show()
    return


def get_log_binning(start, stop, num=100):
    return np.logspace(
        np.log(start),
        np.log(stop),
        num=num,
        endpoint=True,
        base=np.exp(1))


def convert_lambda_to_q(lam, angle):
    angle_conv = np.pi / 180.
    sin_theta_by_2 = np.sin(angle * angle_conv / 2.)
    q = (4. * np.pi / lam) * sin_theta_by_2
    return q


def convert_q_to_lambda(q, angle):
    angle_conv = np.pi / 180.
    sin_theta_by_2 = np.sin(angle * angle_conv / 2.)
    lam = (4. * np.pi / q) * sin_theta_by_2
    return lam


def get_sample_species_info(input_workspace):
    # get sample information: mass, total scattering length, and concentration
    # of each species
    total_stoich = 0.0
    material = mtd[input_workspace].sample().getMaterial().chemicalFormula()
    atom_species = collections.OrderedDict()
    for atom, stoich in zip(material[0], material[1]):
        print(atom.neutron()['tot_scatt_length'])
        b_sqrd_bar = mtd[input_workspace].sample().getMaterial().totalScatterXSection(
        ) / (4. * np.pi)  # <b^2> == sigma_s / 4*pi (in barns)
        atom_species[atom.symbol] = {'mass': atom.mass,
                                     'stoich': stoich,
                                     'b_sqrd_bar': b_sqrd_bar}
        total_stoich += stoich

    for atom, props in atom_species.items(
    ):  # inefficient in py2, but works with py3
        props['concentration'] = props['stoich'] / total_stoich

    return atom_species


def calculate_elastic_self_scattering(input_workspace):
    atom_species = get_sample_species_info(input_workspace)

    # calculate elastic self-scattering term
    elastic_self_term = 0.0
    for species, props in atom_species.items(
    ):  # inefficient in py2, but works with py3
        elastic_self_term += props['concentration'] * props['b_sqrd_bar']

    return elastic_self_term


def calculate_placzek_self_scattering(incident_workspace, output_workspace, L1,
                                      L2, polar, azimuthal=None, detector=None,
                                      parent_workspace=None):
    # constants and conversions
    factor = 1. / scipy.constants.physical_constants['atomic mass unit-kilogram relationship'][0]
    neutron_mass = factor * scipy.constants.m_n

    # get sample information: mass, total scattering length, and concentration
    # of each species
    atom_species = get_sample_species_info(incident_workspace)

    # calculate summation term w/ neutron mass over molecular mass ratio
    summation_term = 0.0
    for species, props in atom_species.items(
    ):  # inefficient in py2, but works with py3
        summation_term += props['concentration'] * \
                          props['b_sqrd_bar'] * neutron_mass / props['mass']

    # get incident spectrum and 1st derivative
    incident_index = 0
    incident_prime_index = 1

    x_lambda = mtd[incident_workspace].readX(incident_index)
    incident = mtd[incident_workspace].readY(incident_index)
    incident_prime = mtd[incident_workspace].readY(incident_prime_index)

    phi_1 = x_lambda * incident_prime / incident

    # Set default Detector Law
    if detector is None:
        detector = {'Alpha': None,
                    'LambdaD': 1.44,
                    'Law': '1/v'}

    # Set detector exponential coefficient alpha
    if detector['Alpha'] is None:
        detector['Alpha'] = 2. * np.pi / detector['LambdaD']

    # Detector law to get eps_1(lambda)
    if detector['Law'] == '1/v':
        c = -detector['Alpha'] / (2. * np.pi)
        x = x_lambda
        detector_law_term = c * x * np.exp(c * x) / (1. - np.exp(c * x))

    eps_1 = detector_law_term

    # Set default azimuthal angle
    if azimuthal is None:
        azimuthal = np.zeros(len(polar))

    # Placzek
    '''
    Original Placzek inelastic correction Ref (for constant wavelength, reactor source):
        Placzek, Phys. Rev v86, (1952), pp. 377-388
    First Placzek correction for time-of-flight, pulsed source (also shows reactor eqs.):
        Powles, Mol. Phys., v6 (1973), pp.1325-1350
    Nomenclature and calculation for this program follows Ref:
         Howe, McGreevy, and Howells, J. Phys.: Condens. Matter v1, (1989), pp. 3433-3451

    NOTE: Powles's Equation for inelastic self-scattering is equal to Howe's Equation for P(theta)
    by adding the elastic self-scattering
    '''
    x_lambdas = np.array([])
    placzek_correction = np.array([])
    for bank, (l2, theta, phi) in enumerate(zip(L2, polar, azimuthal)):
        # variables
        l_total = L1 + l2
        f = L1 / l_total

        angle_conv = np.pi / 180.
        sin_theta = np.sin(theta * angle_conv)
        sin_theta_by_2 = np.sin(theta * angle_conv / 2.)

        term1 = (f - 1.) * phi_1
        term2 = f * eps_1
        term3 = f - 3.

        # per_bank_q = ConvertLambdaToQ(x_lambda,theta)
        inelastic_placzek_self_correction = 2. * \
                                            (
                                                        term1 - term2 + term3) * sin_theta_by_2 * sin_theta_by_2 * summation_term  # See Eq. (A1.14) of
        x_lambdas = np.append(x_lambdas, x_lambda)
        placzek_correction = np.append(
            placzek_correction,
            inelastic_placzek_self_correction)

    if parent_workspace:
        CreateWorkspace(
            DataX=x_lambdas,
            DataY=placzek_correction,
            OutputWorkspace=output_workspace,
            UnitX='Wavelength',
            NSpec=len(polar),
            ParentWorkspace=parent_workspace,
            Distribution=True)
    else:
        CreateWorkspace(
            DataX=x_lambdas,
            DataY=placzek_correction,
            OutputWorkspace=output_workspace,
            UnitX='Wavelength',
            NSpec=len(polar),
            Distribution=True)
    print("Placzek YUnit:", mtd[output_workspace].YUnit())
    print("Placzek distribution:", mtd[output_workspace].isDistribution())

    return mtd[output_workspace]


# -----------------------------------------------------------------------------------------#
# Start Placzek calculations

if '__main__' == __name__:
    # -----------------------------------------------------------------------------------------#
    # Get input parameters
    configfile = sys.argv[1]
    with open(configfile) as handle:
        config = json.loads(handle.read())

    # Get sample and correction info
    sample = config['Sample']
    opts = sample['InelasticCorrection']

    # -----------------------------------------------------------------------------------------#
    # Get incident spectrum
    runs = sample["Runs"].split(',')
    runs = ["%s_%s" % (config["Instrument"], run) for run in runs]
    print("Processing Scan: ", runs[0])

    incident_ws = 'incident_ws'
    get_incident_spectrum_from_monitor(runs[0],
                                       output_workspace=incident_ws)

    # -----------------------------------------------------------------------------------------#
    # Fit incident spectrum
    incident_fit = 'incident_fit'
    fit_type = opts['FitSpectrumWith']
    fit_incident_spectrum(input_workspace=incident_ws,
                          output_workspace=incident_fit,
                          fit_spectrum_with=fit_type,
                          binning_for_fit=opts['LambdaBinningForFit'],
                          binning_for_calc=opts['LambdaBinningForCalc'],
                          plot_diagnostics=opts['PlotFittingDiagnostics'])

    # Set sample info
    SetSampleMaterial(incident_fit, ChemicalFormula=str(sample['Material']))
    calculate_elastic_self_scattering(input_workspace=incident_fit)
    atom_species = get_sample_species_info(incident_fit)

    # Parameters for NOMAD detectors by bank
    L1 = 19.5
    banks = collections.OrderedDict()
    banks[0] = {'L2': 2.01, 'theta': 15.10}
    banks[1] = {'L2': 1.68, 'theta': 31.00}
    banks[2] = {'L2': 1.14, 'theta': 65.00}
    banks[3] = {'L2': 1.11, 'theta': 120.40}
    banks[4] = {'L2': 0.79, 'theta': 150.10}
    banks[5] = {'L2': 2.06, 'theta': 8.60}

    L2 = [x['L2'] for bank, x in banks.iteritems()]
    Polar = [x['theta'] for bank, x in banks.iteritems()]

    parent = 'parent_ws'
    placzek = 'placzek_out'
    Load(Filename=runs[0], OutputWorkspace=parent)
    calculate_placzek_self_scattering(incident_workspace=incident_fit,
                                      output_workspace=placzek,
                                      L1=19.5,
                                      L2=L2,
                                      polar=Polar,
                                      parent_workspace=parent)

    # print(mtd[parent].getNumberHistograms())
    # print(mtd[placzek].getNumberHistograms())
    '''
    ConvertUnits(InputWorkspace=placzek,
             OutputWorkspace=placzek, 
             Target='MomentumTransfer',
             EMode='Elastic')
    '''

    plot = True
    if plot:
        import matplotlib.pyplot as plt

        bank_colors = ['k', 'r', 'b', 'g', 'y', 'c']
        nbanks = range(mtd[placzek].getNumberHistograms())
        for bank, theta in zip(nbanks, Polar):
            x_lambda = mtd[placzek].readX(bank)
            q = convert_lambda_to_q(x_lambda, theta)
            per_bank_placzek = mtd[placzek].readY(bank)
            label = 'Bank: %d at Theta %d' % (bank, int(theta))
            plt.plot(
                q,
                1. +
                per_bank_placzek,
                bank_colors[bank] +
                '-',
                label=label)

        material = ' '.join([symbol +
                             str(int(props['stoich'])) +
                             ' ' for symbol, props in atom_species.iteritems()])
        plt.title('Placzek vs. Q for ' + material)
        plt.xlabel('Q (Angstroms^-1')
        plt.ylabel('1 - P(Q)')
        axes = plt.gca()
        # axes.set_ylim([0.96, 1.0])
        plt.legend()
        plt.show()
