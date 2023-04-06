"""Simulate particle accelerator

Setup and run Monte Carlo simulations of a particle accelerator using the
center of momentum frame for a single process. The simulation procedure
is as follows,

    Use Monte Carlo integration and Monte Carlo sampling to simulate
    scattering events at a fixed energy. The total cross section σ is
    estimated using Monte Carlo integration. The expected number of
    scattering events, N_e, is given by the equation,

    N_e = σ 𝓛_int

    Where σ is the total cross section and 𝓛_int is the time integrated
    luminosity. Asumming the scattering events occur randomly with a
    constant frequency, then the number of events expected within a
    fixed time interval is described by a Poisson distribution. Generate
    a random number N from a Poisson distribution, with mean N_e, to use
    as the number of scattering events during a run. Use rejection
    sampling to generate N samples of random angles from the distribution
    given by the differential cross section. The magnitude of the
    3-momentum of the outgoing particle is,

    p = sqrt(E² - m²)

    where E is the energy of the particle and m is the mass of the
    particle. Data entries are saved as (E, p, θ, φ). The processes being
    simulated are planer (since there are two outgoing particles) and the
    distribution is uniform in the azimuthal angle, so with out loss of
    generality φ = 0.

Note, internal calculations use natural units with ħ=c=1. Energies are
in units of MeV if not otherwise specified.

For more on beam luminosity see Griffiths [1]_ chapter 6 section 1.2 Cross
Sections, the PDG book [2]_ section 31.1 Energy and Luminosity or Thomson
[3]_ chapter 1 section 4 Measurements at particle accelerators.

.. [1] D. J. Griffiths, *Introduction to elementary particles*, 2nd,
   rev. ed. Weinheim: Wiley-VCH, 2008.

.. [2] R.L. Workman et al. (Particle Data Group), Prog. Theor. Exp. Phys.
   2022, 083C01 (2022)

.. [3] M. Thomson, *Modern particle physics*. Cambridge: Cambridge
   University Press, 2013.

"""

import os
import numpy

from . import sampling
from . import utils

class accelerator:
    """Simulate particle accelerator

    Simulate scattering events at a single beam energy or for a seqience
    of beam energies.

    Attributes
    ----------
    L_int : float
        The integrated luminosity of simulation runs.
    dsigma_dOmega : callable
        The differential cross of the interaction process.
    M : float
        The mass of the outgoing particle.

    Methods
    -------
    run(E, fname)
        Run Monte Carlo integration at specified energy.

    run_sequence(E, fname)
        Run Monte Carlo integration for sequence of energies.
    """
    def __init__(self, L_int, dsigma_dOmega, M):
        """Initalize an accelerator

        Setup an accelerator for a specified interaction with a given
        integrated luminosity to be used for simulations.

        Parameters
        ----------
        L_int : float
            The integrated luminosity of simulation runs.
        dsigma_dOmega : callable
            The differential cross of the interaction process.
        M : float
            The mass of the outgoing particle.
        """
        self.L_int = L_int
        self.dsigma_dOmega = dsigma_dOmega
        self.M = M

        # number of samples when using Monte Carlo integration to estimate
        # the total cross section and the maximum of the differential cross
        # section.
        self.samples_MC_int = 1000000


    def _run(self, E, dsigma_max, N_events, fname, append=False):
        """Sample N events from a distribution

        Use Monte Carlo sampling of the differential cross section to
        sample a specified number of events at a fixed energy.

        Parameters
        ----------
        E : float
            The beam energy at wich the simulation is running.
        dsigma_max : float
            Maximum of the distribution given by the differential cross
            section.
        N_events : integer
            Number of events to sample from the distribution.
        fname : string
            Name of data file.
        append : bool
            If True, append output to existing file.
        """
        # generate events
        theta = sampling.monte_carlo_sampling(self.dsigma_dOmega, E, dsigma_max, N_events)

        p_mag = numpy.sqrt(E**2 - self.M**2)

        if append:
            mode = "a"
        else:
            mode = "w"

        os.makedirs(os.path.dirname(fname), exist_ok=True)
        with open(fname, mode) as fout:
            if not append:
                fout.write("# Integrated Luminosity in MeV^2\n")
                fout.write(f"{self.L_int:.7E}\n")
                fout.write("#E, P, theta, phi\n")

            for t in theta:
                fout.write(f"{E:.7E},{p_mag:.7E},{t:.7E},{0:.7E}\n")


    def run(self, E, fname):
        """Run simulation at fixed energy

        Simulate scattering events at a given energy. For more detail
        see the module docstring.

        Parameters
        ----------
        E : float
            The beam energy at wich the simulation is running.
        fname : string
            Name of data file.
        """
        sigma_total, dsigma_max = sampling.monte_carlo_integration(self.dsigma_dOmega, E, self.samples_MC_int)

        # expected number of events
        expected_N_events = sigma_total*self.L_int
        N_events = sampling.get_random_counts(expected_N_events)

        self._run(E, dsigma_max, N_events, fname)


    def run_sequence(self, E, fname):
        """Run simulation for a sequence of energies

        Simulate scattering events for each energy in the sequence and
        split the total integrated luminosity evenly across all of the
        runs. For more detail see the module docstring.

        Parameters
        ----------
        E : array
            Sequence of beam energies at which to run the simulation.
        fname : string
            Name of data file.
        """
        # Luminosity per run
        L_int = self.L_int/len(E)

        first = True
        for i, run_E in enumerate(E):
            utils.progress_bar("Run sequence", i, len(E))

            # estimate total differential cross section
            sigma_total, dsigma_max = sampling.monte_carlo_integration_array(self.dsigma_dOmega, run_E, self.samples_MC_int)

            # expected number of events
            expected_N_events = sigma_total*L_int
            N_events = sampling.get_random_counts(expected_N_events)

            if first:
                first = False
                self._run(run_E, dsigma_max, N_events, fname, append=False)
            else:
                self._run(run_E, dsigma_max, N_events, fname, append=True)
        print()
