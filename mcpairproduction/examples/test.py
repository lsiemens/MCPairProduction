"""Estimat total cross section of muon production

Use totcross to estimate the total interaction cross section of the
neutral weak interaction e⁻ + e⁺ → 𝜇⁻ + 𝜇⁺.

"""

from matplotlib import pyplot
import numpy

from .. import constants
from .. import scattering
from .. import sampling

M_z = constants.M_z

E = numpy.linspace(1, M_z, 3)

def set_fermion(fermion):
    def fun(E, samples):
        theta = samples[:,:,0]
        phi = samples[:,:,1]
        scattering.dsigma_dOmega(E, theta, fermion)
    return fun

data = sampling.monte_carlo_sampling(set_fermion("mu"), E, 10)
print(data)

pyplot.plot(E, scattering.dsigma_dOmega(E, 0, "mu"))
pyplot.show()

