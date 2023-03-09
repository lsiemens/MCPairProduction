"""Validate new code

TODO fill out
Use totcross to estimate the total interaction cross section of the
neutral weak interaction e⁻ + e⁺ → 𝜇⁻ + 𝜇⁺.

"""

from matplotlib import pyplot
import numpy

from .. import accelerator
from .. import scattering
from .. import constants

def function(E, samples):
    return scattering.dsigma_dOmega(E, samples, "mu")

L_int = 21.0E9 # in mBarn^-1
L_int = L_int*constants.hbarc2

LEP = accelerator.accelerator(L_int, function, constants.M_mu)
LEP.run(45E3, "./scratch/run.dat")
