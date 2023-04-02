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

Run = analysis.load(os.path.join(path, E_range_file))

analysis.plot_energy_distribution(Run, n_bins, resolution, Nsigma=Nsigma)
analysis.plot_energy_distribution(Run, n_bins, resolution, theta_range=[0, numpy.pi/2], Nsigma=Nsigma)
analysis.plot_energy_distribution(Run, n_bins, resolution, theta_range=[numpy.pi/2, numpy.pi], Nsigma=Nsigma)

# ----------- FB asymmetry

analysis.plot_FB_asymmetry(Run, n_bins, resolution, Nsigma)
