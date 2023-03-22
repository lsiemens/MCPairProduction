"""Validate new code

TODO fill out
Use totcross to estimate the total interaction cross section of the
neutral weak interaction e⁻ + e⁺ → 𝜇⁻ + 𝜇⁺.

"""

from matplotlib import pyplot
import numpy

from .. import totcross

from .. import sampling
from .. import scattering


# Energy comparison muon
M_z = totcross.M_z

logM_z = numpy.log10(M_z/2)
E_range = (85E3/2, 95E3/2)
samples = 100000
fermion_name = "mu"

def dOmega(samples):
    return numpy.sin(samples)

def function(E, samples):
    A_Z2 = scattering.get_A_Z2(fermion_name)
    return scattering.dsigma_dOmega(E, samples, A_Z2)

E = numpy.linspace(*E_range, 100)
totcross_cross, _ = totcross.sigma_estimate(E, fermion_name, samples)
sampling_cross, _ = sampling.monte_carlo_integration_array(function, E, dOmega, samples)
analytic = totcross.sigma_analytic(E, fermion_name)
pyplot.plot(E, totcross_cross, label="MC integral from totcross")
pyplot.plot(E, sampling_cross, label="MC integral from sampling")
pyplot.plot(E, analytic, label="MC analytic")
pyplot.title("Comparison")
pyplot.legend()
pyplot.show()


# angular distribution comparison muon
E = 45E3
samples = 10000

theta = numpy.linspace(0, numpy.pi, 1000)

_, maximum = sampling.monte_carlo_integration(function, E, dOmega, samples)
sampling_diffcross = sampling.monte_carlo_sampling(function, E, maximum, samples)

totcross_dist = totcross.dsigma_dOmega(E, theta, fermion_name)
pyplot.hist(sampling_diffcross, bins=50, alpha=0.5, density=True)
pyplot.plot(theta, totcross_dist/numpy.mean(totcross_dist*numpy.pi))
pyplot.show()
