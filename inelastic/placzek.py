import sys
import json
from mantid import mtd
from mantid.simpleapi import *

from inelastic.incident_spectrum import GetIncidentSpectrumFromMonitor, FitIncidentSpectrum


#-------------------------------------------------------------------------
# Placzek - 1st order inelastic correction

def plotPlaczek(x, y, fit, fit_prime, title=None):
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


def GetLogBinning(start, stop, num=100):
    return np.logspace(
        np.log(start),
        np.log(stop),
        num=num,
        endpoint=True,
        base=np.exp(1))


def ConvertLambdaToQ(lam, angle):
    angle_conv = np.pi / 180.
    sin_theta_by_2 = np.sin(angle * angle_conv / 2.)
    q = (4. * np.pi / lam) * sin_theta_by_2
    return q


def ConvertQToLambda(q, angle):
    angle_conv = np.pi / 180.
    sin_theta_by_2 = np.sin(angle * angle_conv / 2.)
    lam = (4. * np.pi / q) * sin_theta_by_2
    return lam


def CalculateElasticSelfScattering(InputWorkspace):
    # get sample information: mass, total scattering length, and concentration
    # of each species
    total_stoich = 0.0
    material = mtd[InputWorkspace].sample().getMaterial().chemicalFormula()
    atom_species = collections.OrderedDict()
    for atom, stoich in zip(material[0], material[1]):
        print(atom.neutron()['tot_scatt_length'])
        b_sqrd_bar = mtd[InputWorkspace].sample().getMaterial().totalScatterXSection(
        ) / (4. * np.pi)  # <b^2> == sigma_s / 4*pi (in barns)
        atom_species[atom.symbol] = {'mass': atom.mass,
                                     'stoich': stoich,
                                     'b_sqrd_bar': b_sqrd_bar}
        total_stoich += stoich

    for atom, props in atom_species.items(
    ):  # inefficient in py2, but works with py3
        props['concentration'] = props['stoich'] / total_stoich

    # calculate elastic self-scattering term
    elastic_self_term = 0.0
    for species, props in atom_species.items(
    ):  # inefficient in py2, but works with py3
        elastic_self_term += props['concentration'] * props['b_sqrd_bar']

    return elastic_self_term


def CalculatePlaczekSelfScattering(
        IncidentWorkspace,
        ParentWorkspace,
        OutputWorkspace,
        L1,
        L2,
        Polar,
        Azimuthal=None,
        Detector=None):

    # constants and conversions
    neutron_mass = m_n / \
        physical_constants['atomic mass unit-kilogram relationship'][0]

    # get sample information: mass, total scattering length, and concentration
    # of each species
    total_stoich = 0.0
    material = mtd[IncidentWorkspace].sample().getMaterial().chemicalFormula()
    atom_species = collections.OrderedDict()
    for atom, stoich in zip(material[0], material[1]):
        print(atom.neutron()['tot_scatt_length'])
        b_sqrd_bar = mtd[IncidentWorkspace].sample().getMaterial(
        ).totalScatterXSection() / (4. * np.pi)  # <b^2> == sigma_s / 4*pi (in barns)
        atom_species[atom.symbol] = {'mass': atom.mass,
                                     'stoich': stoich,
                                     'b_sqrd_bar': b_sqrd_bar}
        total_stoich += stoich

    for atom, props in atom_species.items(
    ):  # inefficient in py2, but works with py3
        props['concentration'] = props['stoich'] / total_stoich

    # calculate summation term w/ neutron mass over molecular mass ratio
    summation_term = 0.0
    for species, props in atom_species.items(
    ):  # inefficient in py2, but works with py3
        summation_term += props['concentration'] * \
            props['b_sqrd_bar'] * neutron_mass / props['mass']

    # get incident spectrum and 1st derivative
    incident_index = 0
    incident_prime_index = 1

    x_lambda = mtd[IncidentWorkspace].readX(incident_index)
    incident = mtd[IncidentWorkspace].readY(incident_index)
    incident_prime = mtd[IncidentWorkspace].readY(incident_prime_index)

    phi_1 = x_lambda * incident_prime / incident

    # Set default Detector Law
    if Detector is None:
        Detector = {'Alpha': None,
                    'LambdaD': 1.44,
                    'Law': '1/v'}

    # Set detector exponential coefficient alpha
    if Detector['Alpha'] is None:
        Detector['Alpha'] = 2. * np.pi / Detector['LambdaD']

    # Detector law to get eps_1(lambda)
    if Detector['Law'] == '1/v':
        c = -Detector['Alpha'] / (2. * np.pi)
        x = x_lambda
        detector_law_term = c * x * np.exp(c * x) / (1. - np.exp(c * x))

    eps_1 = detector_law_term

    # Set default azimuthal angle
    if Azimuthal is None:
        Azimuthal = np.zeros(len(Polar))

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
    for bank, (l2, theta, phi) in enumerate(zip(L2, Polar, Azimuthal)):
        # variables
        L_total = L1 + l2
        f = L1 / L_total

        angle_conv = np.pi / 180.
        sin_theta = np.sin(theta * angle_conv)
        sin_theta_by_2 = np.sin(theta * angle_conv / 2.)

        term1 = (f - 1.) * phi_1
        term2 = f * eps_1
        term3 = f - 3.

        #per_bank_q = ConvertLambdaToQ(x_lambda,theta)
        inelastic_placzek_self_correction = 2. * \
            (term1 - term2 + term3) * sin_theta_by_2 * sin_theta_by_2 * summation_term  # See Eq. (A1.14) of
        x_lambdas = np.append(x_lambdas, x_lambda)
        placzek_correction = np.append(
            placzek_correction,
            inelastic_placzek_self_correction)

    CreateWorkspace(
        DataX=x_lambdas,
        DataY=placzek_correction,
        OutputWorkspace=OutputWorkspace,
        UnitX='Wavelength',
        NSpec=len(Polar),
        ParentWorkspace=ParentWorkspace,
        Distribution=True)
    print("Placzek YUnit:", mtd[OutputWorkspace].YUnit())
    print("Placzek distribution:", mtd[OutputWorkspace].isDistribution())

    return mtd[OutputWorkspace]

#-----------------------------------------------------------------------------------------#
# Start Placzek calculations

if '__main__' == __name__:
    #-----------------------------------------------------------------------------------------#
    # Get input parameters
    configfile = sys.argv[1]
    with open(configfile) as handle:
        config = json.loads(handle.read())

    # Get sample and correction info
    sample = config['Sample']
    opts = sample['InelasticCorrection']


    #-----------------------------------------------------------------------------------------#
    # Get incident spectrum
    print("Processing Scan: ", sample['Runs'][0])

    incident_ws = 'incident_ws'
    lam_binning = opts['LambdaBinningForFit']
    GetIncidentSpectrumFromMonitor(sample['Runs'][0],
                                   OutputWorkspace=incident_ws,
                                   LambdaBinning=lam_binning)

    #-----------------------------------------------------------------------------------------#
    # Fit incident spectrum
    incident_fit = 'incident_fit'
    fit_type = opts['FitSpectrumWith']
    FitIncidentSpectrum(InputWorkspace=incident_ws,
                        OutputWorkspace=incident_fit,
                        FitSpectrumWith=fit_type,
                        BinningForFit=opts['LambdaBinningForFit'],
                        BinningForCalc=opts['LambdaBinningForCalc'],
                        PlotDiagnostics=opts['PlotFittingDiagnostics'])

    # Set sample info
    SetSampleMaterial(incident_fit, ChemicalFormula=str(sample['Material']))
    atom_species = GetSamplePropsForInelasticCorr(incident_fit)

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

    CalculatePlaczekSelfScattering(IncidentWorkspace=incident_fit,
                                   OutputWorkspace='placzek_out',
                                   L1=19.5,
                                   L2=L2,
                                   Polar=Polar)

    if plot:
        import matplotlib.pyplot as plt
        bank_colors = ['k', 'r', 'b', 'g', 'y', 'c']
        nbanks = range(mtd['placzek_out'].getNumberHistograms())
        for bank, theta in zip(nbanks, Polar):
            q = mtd['placzek_out'].readX(bank)
            per_bank_placzek = mtd['placzek_out'].readY(bank)
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
        axes.set_ylim([0.96, 1.0])
        plt.legend()
        plt.show()
