#!/usr/bin/env python3
from matplotlib import pyplot
import numpy
import os

from mcpairproduction import scattering
from mcpairproduction import analysis

path = "./scratch/"

resolution = 1000
n_bins = 32
Nsigma = 2 # Size of error bars in sigma

with open(os.path.join(path, "data_index"), "r") as fin:
    E_files = fin.read()
    E_files = E_files.split("\n")

E_files, E_range_file = E_files[:-1], E_files[-1]

Runs_fixed_E = [analysis.load(os.path.join(path, file)) for file in E_files]

# ----------- Observations of angular distribution at fixed energy
# differential cross section functions

for Run in Runs_fixed_E:
    analysis.plot_angular_distribution(Run, n_bins, resolution, Nsigma)

# -------------- E range run
# total cross section functions

def sigma_analytic(E):
    A_tot2 = scattering.get_A_tot2_integrated()
    return scattering.sigma_analytic(E, A_tot2)

Run = analysis.load(os.path.join(path, E_range_file))
E_data, p_data, theta_data, phi_data, L_int = Run

E_range = (numpy.min(E_data), numpy.max(E_data))
E_theory = numpy.linspace(*E_range, resolution)

sigma_total = sigma_analytic(E_theory)

events_per_bin = len(E_data)/n_bins
#normalization = events_per_bin/numpy.mean(sigma_total)
normalization = events_per_bin/numpy.mean(sigma_analytic(numpy.array(list(set(E_data)))))

hist, bin_edges = numpy.histogram(E_data, bins=n_bins, range=E_range)
hist_err = Nsigma*numpy.sqrt(hist)

bin_width = bin_edges[1] - bin_edges[0]
pyplot.bar(bin_edges[:-1], hist, bin_width, align="edge", yerr=hist_err)

pyplot.plot(E_theory, normalization*sigma_total, "k--")
pyplot.show()

# ----------- FB asymmetry
# integrated cross section functions

def sigma_forward(E):
    A_2_forward = scattering.get_A_tot2_integrated(theta_range=(0, numpy.pi/2))
    return scattering.sigma_analytic(E, A_2_forward)

def sigma_back(E):
    A_2_back = scattering.get_A_tot2_integrated(theta_range=(numpy.pi/2, numpy.pi))
    return scattering.sigma_analytic(E, A_2_back)

front_mask = theta_data <= numpy.pi/2

Fhist, bin_edges = numpy.histogram(E_data[front_mask], bins=n_bins, range=E_range)
Fhist_err = Nsigma*numpy.sqrt(Fhist)

FB_asymmetry = 2*Fhist/hist - 1 # (F - B) / T, F + B = T

asymmetry_err = (2*Fhist/hist)*numpy.sqrt((Fhist_err/Fhist)**2 + (hist_err/hist)**2)
FB_asymmetry_analytic = 2*sigma_forward(E_theory)/sigma_analytic(E_theory) - 1

bin_width = bin_edges[1] - bin_edges[0]
pyplot.bar(bin_edges[:-1], FB_asymmetry, bin_width, align="edge", yerr=asymmetry_err)

pyplot.plot(E_theory, FB_asymmetry_analytic, "k--")
pyplot.title("FB_asymmetry")
pyplot.show()

pyplot.bar(bin_edges[:-1], Fhist, bin_width, align="edge", yerr=Fhist_err)

pyplot.plot(E_theory, normalization*sigma_forward(E_theory), "k--")
pyplot.title("sigma forward")
pyplot.show()

pyplot.bar(bin_edges[:-1], hist - Fhist, bin_width, align="edge", yerr=Nsigma*numpy.sqrt(hist - Fhist))

pyplot.plot(E_theory, normalization*sigma_back(E_theory), "r-.")
pyplot.title("sigma back")
pyplot.show()
