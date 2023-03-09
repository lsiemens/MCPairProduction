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
    return scattering.dsigma_dOmega(E, samples, fermion_name)

E = numpy.linspace(*E_range, 100)
totcross_cross, _ = totcross.sigma_estimate(E, fermion_name, samples)
sampling_cross, _ = sampling.monte_carlo_integration(function, E, dOmega, samples)
analytic = totcross.sigma_analytic(E, fermion_name)
pyplot.plot(E, totcross_cross, label="MC integral from totcross")
pyplot.plot(E, sampling_cross, label="MC integral from sampling")
pyplot.plot(E, analytic, label="MC analytic")
pyplot.title("Comparison")
pyplot.legend()
pyplot.show()


# angular distribution comparison muon
num_E = 3
E_range = (1E3, 90E3)
E = numpy.linspace(*E_range, num_E)
samples = 10000

theta = numpy.linspace(0, numpy.pi, 1000)

_, maximum = sampling.monte_carlo_integration(function, E, dOmega, samples)
sampling_diffcross = sampling.monte_carlo_sampling(function, E, maximum, samples)

for i in range(num_E):
    totcross_dist = totcross.dsigma_dOmega(E[i], theta, fermion_name)
    pyplot.hist(sampling_diffcross[:, i], bins=50, alpha=0.3, density=True)
    pyplot.plot(theta, totcross_dist/numpy.mean(totcross_dist*numpy.pi))
pyplot.show()
