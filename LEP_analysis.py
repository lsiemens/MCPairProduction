#!/usr/bin/env python3
from matplotlib import pyplot
import numpy
import os

from mcpairproduction import scattering
from mcpairproduction import analysis

path = "./scratch/"

resolution = 1000
n_bins = 25
Nsigma = 2 # errors

E_files = ["RUN_42GeV.dat", "RUN_45GeV.dat", "RUN_48GeV.dat"]
E_range_file = "RUN_E_range.dat"

Runs_fixed_E = [analysis.load(os.path.join(path, file)) for file in E_files]
Run_E_range = analysis.load(os.path.join(path, E_range_file))

# ----------- Observations of angular distribution at fixed energy
Run = Runs_fixed_E[0]
E_data, p_data, theta_data, phi_data, L_int = Run

E_data = E_data[0]
theta_theory = numpy.linspace(0, numpy.pi, resolution)

# differential cross section functions

def dsigma_dOmega(E, theta):
    A_tot2 = scattering.get_A_tot2("mu")
    return scattering.dsigma_dOmega(E, theta, A_tot2)

def dsigma_dOmega_gamma(E, theta):
    A2 = scattering.get_A_gamma2("mu")
    return scattering.dsigma_dOmega(E, theta, A2)

def dsigma_dOmega_Z(E, theta):
    A2 = scattering.get_A_Z2("mu")
    return scattering.dsigma_dOmega(E, theta, A2)

def dsigma_dOmega_interference(E, theta):
    A2 = scattering.get_A_cross2("mu")
    return scattering.dsigma_dOmega(E, theta, A2)

# calibrating the differential cross section to compair with the data
ds_dO_total = dsigma_dOmega(E_data, theta_theory)
ds_dO_gamma = dsigma_dOmega_gamma(E_data, theta_theory)
ds_dO_Z = dsigma_dOmega_Z(E_data, theta_theory)
ds_dO_interference = dsigma_dOmega_interference(E_data, theta_theory)

events_per_bin = len(theta_data)/n_bins
normalization = events_per_bin/numpy.mean(ds_dO_total)
#ds_dO_events = events_per_bin*ds_dO_total/numpy.mean(ds_dO_total)

# make histogram of angular distribution
hist, bin_edges = numpy.histogram(theta_data, bins=n_bins, range=(0, numpy.pi))
hist_err = Nsigma*numpy.sqrt(hist) # 2 sigma error from poisson statistics

bin_width = bin_edges[1] - bin_edges[0]
pyplot.bar(bin_edges[:-1], hist, bin_width, align="edge", yerr=hist_err)

pyplot.plot(theta_theory, normalization*ds_dO_total)
pyplot.plot(theta_theory, normalization*ds_dO_gamma)
pyplot.plot(theta_theory, normalization*ds_dO_Z)
pyplot.plot(theta_theory, normalization*ds_dO_interference)
pyplot.show()
