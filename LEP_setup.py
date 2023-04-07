#!/usr/bin/env python3
"""Setup LEP simulation

This script simulates particle scattering by the process
e⁻ + e⁺ →  𝜇⁻ + 𝜇⁺ at LEP. For the integrated luminosity and beam energy
of the LEP see [1]_.

.. [1] Arduini, G., R. Assmann, R. Bailey, A. Butterworth, P. Collier,
       K. Cornelis, S. Fartoukh, et al. "Electron-Positron Collisions at
       209 GeV in LEP." In PACS2001. Proceedings of the 2001 Particle
       Accelerator Conference (Cat. No.01CH37268), 1:356-358 vol.1. IEEE,
       2001.

"""

from matplotlib import pyplot
import numpy
import os

from mcpairproduction import accelerator
from mcpairproduction import scattering
from mcpairproduction import constants

def dsigma_dOmega(E, theta):
    """Differential scattering cross section

    Parameters
    ----------
    E : array
        The beam energy.
    theta : array
        The angle of the outgoing muon.
    """
    A_tot2 = scattering.get_A_tot2("mu")
    return scattering.dsigma_dOmega(E, theta, A_tot2)

path = "./scratch/"

# format (Max beam Energy in MeV, Integrated luminosity in mbarn^-1)
LEP_operation = {1994:(4.56E4,  6.4E10),
                 1995:(7.00E4,  4.7E10),
                 1996:(8.60E4,  2.5E10),
                 1997:(9.20E4,  7.5E10),
                 1998:(9.45E4,  2.00E11),
                 1999:(1.010E5, 2.54E11),
                 2000:(1.045E5, 2.33E11)}

year = 2000
Max_E, L_int = LEP_operation[year]
L_int = L_int*constants.hbarc2 # in MeV^2

LEP = accelerator.accelerator(L_int, dsigma_dOmega, constants.M_mu)

# Set the energy for fixed energy runs
E_values = [2.6E4, 4.0E4, 5.1E4]

fnames = [f"RUN_{int(E/1E3)}GeV.dat" for E in E_values]
for E, fname in zip(E_values, fnames):
    fname = f"RUN_{int(E/1E3)}GeV.dat"
    LEP.run(E, os.path.join(path, fname))

# The number of Energy steps should be much larger than the number of
# bins in the analysis to avoid bias.
N_steps = 1024
E_range = (10*constants.M_mu, Max_E)
E = numpy.linspace(*E_range, N_steps)

fname = "RUN_E_range.dat"
LEP.run_sequence(E, os.path.join(path, fname))

# record list of data files produced for use by the analysis script
fnames += [fname]
with open(os.path.join(path, "data_index"), "w") as fout:
    content = "\n".join(fnames)
    fout.write(content)
