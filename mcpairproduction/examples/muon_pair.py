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
from .. import analysis

def function(E, samples):
    A_tot2 = scattering.get_A_cross2("mu")
    return scattering.dsigma_dOmega(E, samples, A_tot2)

def analytic_integral(E):
    from .. import totcross
    return totcross.sigma_analytic(E, "mu")

path = "./scratch/run.dat"
n_bins = 25

L_int = 100*21.0E9/12 # in mBarn^-1
L_int = L_int*constants.hbarc2

#E = 45E3
E = constants.M_z/2 - 1E3

LEP = accelerator.accelerator(L_int, function, constants.M_mu)
LEP.run(E, path)

E_, p_mag_, theta, phi_ = analysis.load(path)

hist, bin_edges = numpy.histogram(theta, bins=n_bins)
err = numpy.sqrt(hist)

print("fig1")
theta = numpy.linspace(0, numpy.pi)
dist = function(E, theta)
dist = dist*numpy.mean(hist)/numpy.mean(dist)
pyplot.bar(bin_edges[:-1], hist, width=numpy.diff(bin_edges), yerr=err, align="edge")
pyplot.plot(theta, dist, "r")
pyplot.show()

print("fig2")
symm = (hist - hist[::-1])/2
symm_err = err + err[::-1]
pyplot.bar(bin_edges[:-1], symm, width=numpy.diff(bin_edges), yerr=symm_err, align="edge")
pyplot.show()


print("fig3 E dist")
# E sequence
path2 = "./scratch/run2.dat"
E_range = (1E3, 90E3)
E = numpy.linspace(*E_range, 1000)
LEP.run_sequence(E, path2)

E_, p_mag_, theta, phi_ = analysis.load(path2)

bin_edges = numpy.linspace(*E_range, n_bins + 1)
hist, bin_edges = numpy.histogram(E_, bins=bin_edges)
err = numpy.sqrt(hist)

theta = numpy.linspace(0, numpy.pi)
pyplot.bar(bin_edges[:-1], hist, width=numpy.diff(bin_edges), yerr=err, align="edge")

dist = analytic_integral(E)*LEP.L_int/n_bins
pyplot.plot(E, dist, "r")
pyplot.show()

print("END")
