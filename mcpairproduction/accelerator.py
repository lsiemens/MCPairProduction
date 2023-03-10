""" luminosity / write data / ...

Estimate total cross sections of neutral weak interactions.

Estimates the total interaction cross section for electron positron
scattering mediated by the Z boson e⁻ + e⁺ →  f + f̄. The estimate of
the total cross section is computed from the differential cross section
using Monte Carlo integration.

Note, internal calculations use natural units with ħ=c=1. Energies are
in units of MeV if not otherwise specified.

The differentail cross section is from Griffiths [1]_ Example 9.5
Electron-Positron scattering near the Z pole, with *Review of Particle
Physics* [2]_ as an additional refrence and source for physical constants.
Equations for calculating the weak coupling to the Z boson are from
Thomson [3]_ section 15.3.1

.. [1] D. J. Griffiths, *Introduction to elementary particles*, 2nd,
   rev. ed. Weinheim: Wiley-VCH, 2008.

.. [2] R.L. Workman et al. (Particle Data Group), Prog. Theor. Exp. Phys.
   2022, 083C01 (2022)

.. [3] M. Thomson, *Modern particle physics*. Cambridge: Cambridge
   University Press, 2013.

"""

import numpy

from . import sampling

class accelerator:
    def __init__(self, L_int, dsigma_dOmega, M):
        self.L_int = L_int
        self.dsigma_dOmega = dsigma_dOmega
        self.M = M

        def dOmega(samples):
            return numpy.sin(samples)

        self.dOmega = dOmega

        # number of samples when using Monte Carlo integration to estimate
        # the total cross section and the maximum of the differential cross
        # section.
        self.samples_MC_int = 1000


    def _run(self, E, dsigma_max, N_events, path, append=False):
        # generate events
        theta = sampling.monte_carlo_sampling(self.dsigma_dOmega, E, dsigma_max, N_events)

        p_mag = numpy.sqrt(E**2 - self.M**2)

        if append:
            mode = "a"
        else:
            mode = "w"

        with open(path, mode) as fout:
            if not append:
                fout.write("#E, P, theta, phi\n")

            for t in theta:
                fout.write(f"{E:.5E},{p_mag:.5E},{t:.5E},{0:.5E}\n")


    def run(self, E, path, append=False):
        sigma_total, dsigma_max = sampling.monte_carlo_integration(self.dsigma_dOmega, E, self.dOmega, self.samples_MC_int)

        # expected number of events
        expected_N_events = sigma_total*self.L_int
        N_events = sampling.get_random_counts(expected_N_events)

        self._run(E, dsigma_max, N_events, path, append)


    def run_sequence(self, E, path):
        """Collect data at a sequence of energies
        """
        sigma_total, dsigma_max = sampling.monte_carlo_integration_array(self.dsigma_dOmega, E, self.dOmega, self.samples_MC_int)

        # Luminosity per run
        L_int = self.L_int/len(E)

        # expected number of events
        expected_N_events = sigma_total*L_int
        N_events = sampling.get_random_counts(expected_N_events)

        first = True
        for run_E, run_dsigma_max, run_N_events in zip(E, dsigma_max, N_events):
            if first:
                first = False
                self._run(run_E, run_dsigma_max, run_N_events, path, append=False)
            else:
                self._run(run_E, run_dsigma_max, run_N_events, path, append=True)
