# package imports
from total_scattering.inelastic.incident_spectrum import (
    fitCubicSpline, fitCubicSplineViaMantidSplineSmoothing,
    fitCubicSplineWithGaussConv, FitIncidentSpectrum)
from total_scattering.utils import random_name
from tests import TEST_DATA_DIR

# third-party imports
from mantid.simpleapi import (CreateWorkspace, DeleteWorkspace,
                              DeleteWorkspaces, LoadAscii, mtd, Rebin)
import numpy as np

# standard imports
from pathlib import Path
import unittest


class TestIncidentSpectrum(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        # Load experimental incident spectrum
        cls.incident_spectrum = random_name()
        file_path = str(Path(TEST_DATA_DIR) /
                        'NOM_161959_incident_spectrum.dat')
        LoadAscii(Filename=file_path,
                  OutputWorkspace=cls.incident_spectrum,
                  Unit='Wavelength')
        x = mtd[cls.incident_spectrum].readX(0)
        y = mtd[cls.incident_spectrum].readY(0)
        # fit to 10 degree polynomial. Eeplace spectrum with the polynomial,
        # necessary to precisely evaluate the derivatives yielded
        # by the interpolators
        poly = np.polynomial.Polynomial.fit(x, y, 10)
        CreateWorkspace(x, poly(x), UnitX='Wavelength',
                        OutputWorkspace=cls.incident_spectrum)

        # x and y values used to fit the incident spectrum
        cls.bin_width = 0.05
        cls.BinningForFit = [min(x)+0.05, cls.bin_width, max(x) - 0.1]
        cls.spectrum_to_fit = random_name()
        Rebin(InputWorkspace=cls.incident_spectrum,
              Params=cls.BinningForFit,
              OutputWorkspace=cls.spectrum_to_fit)
        # readX() and readY() return references, thus return a copy
        cls.x_fit = np.array(mtd[cls.spectrum_to_fit].readX(0))
        cls.y_fit = np.array(mtd[cls.spectrum_to_fit].readY(0))

        # x values where to evaluate the fitted incident spectrum
        bin_width = (max(x) - min(x)) / len(x)
        cls.BinningForCalc = [min(x)+0.1, bin_width, max(x) - 0.15]
        temp = Rebin(InputWorkspace=cls.incident_spectrum,
                     Params=cls.BinningForCalc,
                     OutputWorkspace=random_name())
        x = np.array(temp.readX(0))
        cls.y = poly(x)
        cls.y_prime = poly.deriv()(x)
        cls.y_prime2 = poly.deriv(2)(x)
        cls.x = x
        DeleteWorkspace(temp)

    @classmethod
    def tearDownClass(cls) -> None:
        DeleteWorkspaces([cls.incident_spectrum, cls.spectrum_to_fit])

    def assess_fit(self, y, y_prime, y_prime2, correlations):
        # the correlation between the obtained and expected values will
        # degrade with the order of the derivative
        for obtained, expected, corr in zip(
                [y, y_prime, y_prime2],
                [self.y, self.y_prime, self.y_prime2],
                correlations):
            corr_obtained = np.corrcoef(obtained, expected)[0][1]
            self.assertTrue(corr_obtained > corr,
                            f'correlation {corr_obtained} lower than {corr}')

    def test_fitCublicSpline(self):
        fit = fitCubicSpline(self.x_fit, self.y_fit, self.x)
        self.assess_fit(*fit, [0.999, 0.99, 0.97])

    def test_fitCubicSplineViaMantidSplineSmoothing(self):
        fit = fitCubicSplineViaMantidSplineSmoothing(self.spectrum_to_fit,
                                                     self.BinningForCalc,
                                                     MaxNumberOfBreaks=32)
        # the correlation between the obtained values for the incident
        # spectrum and the expected one is very high (0.9999) but the
        # correlations for the first and second derivatives are poor
        self.assess_fit(*fit, [0.999, 0.80, 0.10])

    def test_fitCubicSplineWithGaussConv(self):
        fit = fitCubicSplineWithGaussConv(self.x_fit, self.y_fit, self.x,
                                          sigma=4 * self.bin_width)
        self.assess_fit(*fit, [0.999, 0.99, 0.97])

    def test_FitIncidentSpectrum(self):
        r"""Just ensure that it runs with no exception"""
        output_workspace = random_name()
        for interpolator in ['CubicSpline',
                             'CubicSplineViaMantid',
                             'HowellsFunction',
                             'GaussConvCubicSpline']:
            FitIncidentSpectrum(self.incident_spectrum, output_workspace,
                                FitSpectrumWith=interpolator,
                                BinningForFit=self.BinningForFit,
                                BinningForCalc=self.BinningForCalc)
        DeleteWorkspace(output_workspace)
