#!/usr/bin/env python3

"""
Supplemental figures for the presentation
"""

import sys
sys.path.append("../..")

from matplotlib import pyplot
import matplotlib.patches as patches
import numpy
from mcpairproduction import scattering
from mcpairproduction import sampling

save_fig = True
scale = 1.5
seed = 2
pyplot.rcParams.update({"savefig.bbox":"tight",
                        "savefig.dpi":300})

resolution = 100
E = 20E3 # in MeV

function = lambda theta: scattering.dsigma_dOmega(E, theta, scattering.get_A_tot2())

domain = numpy.linspace(0, numpy.pi, resolution)

data = function(domain)
normalization = 1/numpy.max(data)

#  ------------    Monte Carlo Integraion
Npoints = 20
example_index = 5

samples, _ = sampling.get_random_samples((Npoints,), seed=seed)
sampled_func = function(samples)

sampled_mean = numpy.mean(sampled_func)

marker_size = 4*pyplot.rcParams["lines.markersize"]**2

for i,v in enumerate(samples):
    print(i, v/numpy.pi)

fig, ax = pyplot.subplots()
ax.plot(domain, normalization*data, "k-", label="$f(x)$")
ax.scatter(samples, normalization*sampled_func, marker="2", color="b", s=marker_size, label="Value at sample")
ax.hlines(normalization*sampled_mean, 0, numpy.pi, linestyle="--", color="k", label="$f_{avg}$")
ax.vlines(samples[example_index], 0, normalization*function(samples[example_index]), linestyle=":", color="b", label="Sample $x_i$")
ax.set_ylim(0, None)
ax.set_xlim(0, numpy.pi)
ax.set_xlabel("Domain")
ax.tick_params(axis="both", labeltop=False, labelbottom=False, labelleft=False, labelright=False)
ax.legend()

if save_fig:
    size = fig.get_size_inches()
    size = size/scale
    fig.set_size_inches(size)
    pyplot.savefig("MCIntegral.png")
pyplot.show()


# -------------    Monte Carlo Sampling
Npoints = 100

samples = sampling.get_random_rejection_samples((Npoints,), seed=seed)
theta = samples[0, :]
v = samples[1, :]

mask = (v <= normalization*function(theta))
mask_reject = numpy.logical_not(mask)

fig, ax = pyplot.subplots()
ax.plot(domain, normalization*data, "k-", label="$f(x)/f_{max}$")
ax.scatter(theta[mask], v[mask], marker="o", color="g", label="Accepted")
ax.scatter(theta[mask_reject], v[mask_reject], marker="X", color="r", label="Rejected")
ax.set_ylim(0, None)
ax.set_xlim(0, numpy.pi)
ax.set_xlabel("Domain")
ax.tick_params(axis="both", labeltop=False, labelbottom=False, labelleft=False, labelright=False)
ax.legend(framealpha=1)


if save_fig:
    size = fig.get_size_inches()
    size = size/scale
    fig.set_size_inches(size)
    pyplot.savefig("MCSampling.png")
pyplot.show()
