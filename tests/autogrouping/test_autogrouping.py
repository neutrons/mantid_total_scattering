import unittest
import numpy as np
import random

from total_scattering.autogrouping.autogrouping import gather_fitparameters

from mantid.simpleapi import \
    CreateEmptyTableWorkspace, \
    FitPeaks, \
    GeneratePeaks, \
    mtd
from Calibration.tofpd import diagnostics

DIAMOND_PEAKS = (0.8920, 1.0758, 1.2615)
PARAMETERS = ["Mixing", "Intensity", "PeakCentre", "FWHM"]


class TestAutogrouping(unittest.TestCase):

    def setUp(self):
        return

    def test_gatherfit(self):
        peakpositions = np.asarray(DIAMOND_PEAKS)
        peakwindows = diagnostics.get_peakwindows(peakpositions)

        # Create table workspace to hold peak parameters
        params = CreateEmptyTableWorkspace()
        params.addColumn("int", "spectrum")
        params.addColumn("double", "Mixing")
        params.addColumn("double", "Intensity")
        params.addColumn("double", "PeakCentre")
        params.addColumn("double", "FWHM")
        params.addColumn("double", "A0")
        params.addColumn("double", "A1")
        params.addColumn("double", "chi2")

        # Generate peaks with some random variation
        n = 10
        sigma = 0.3
        for i in range(0, n):
            # randomly choose a diamond peak center
            peak = random.randint(0, 2)

            x0 = random.gauss(peakpositions[peak], sigma)
            intensity = random.gauss(1.0, 0.25)
            fwhm = random.gauss(1.0, 0.25)
            mixing = random.gauss(0.6, 0.2)
            params.addRow([i, mixing, intensity, x0, fwhm, 0.0, 0.0, 0.01])

        genpeaks = GeneratePeaks(PeakParametersWorkspace=params, PeakType='PseudoVoigt',
                                 BackgroundType='Linear (A0, A1)',
                                 BinningParameters='0.1,0.001,10', NumberWidths=2)
        genpeaks.setYUnit("Counts")

        # Fit the peaks
        FitPeaks(InputWorkspace='genpeaks',
                 OutputWorkspace='output',
                 PeakFunction='PseudoVoigt',
                 PeakParameterNames="Mixing",
                 PeakParameterValues="0.6",
                 RawPeakParameters=True,
                 ConstrainPeakPositions=False,
                 PeakCenters=peakpositions,
                 FitWindowBoundaryList=peakwindows,
                 FittedPeaksWorkspace='fitted',
                 OutputPeakParametersWorkspace='parameters',
                 OutputParameterFitErrorsWorkspace='fiterrors')

        # Convert the parameter workspace to array used by clustering algorithm
        result, mask = gather_fitparameters('parameters', PARAMETERS, None,
                                            peakpositions)
        result = result[..., 1:]  # remove wsindex column

        npeaks = len(peakpositions)
        nparams = len(PARAMETERS)
        self.assertEqual(result.shape, (n, npeaks * nparams))

        # Check each peak parameter against fitpeaks result
        for i in range(n):
            for j in range(npeaks):
                row = mtd['parameters'].row(i * npeaks + j)
                for k in range(nparams):
                    self.assertEqual(result[i, nparams * j + k], row[PARAMETERS[k]])


if __name__ == '__main__':
    unittest.main()
