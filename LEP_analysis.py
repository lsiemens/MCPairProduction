#!/usr/bin/env python3
from matplotlib import pyplot
import numpy
import os

from mcpairproduction import scattering
from mcpairproduction import analysis

path = "./scratch/"

fig_path = "./presentation/Figures"
save_fig = True
scale = 1.4
pyplot.rcParams.update({"savefig.bbox":"tight",
                        "savefig.dpi":400})

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
    fig, ax = pyplot.subplots()
    analysis.plot_angular_distribution(ax, Run, n_bins, resolution, Nsigma)

    if save_fig:
        E, _, _, _, _ = Run
        E = E[0]
        size = fig.get_size_inches()
        size = size/scale
        fig.set_size_inches(size)
        pyplot.savefig(os.path.join(fig_path, f"{E/1E3:.0f}GeV_angular_distribution.png"))
    pyplot.show()

# -------------- E range run

Run = analysis.load(os.path.join(path, E_range_file))
comparisons = [("tot", "k-", False), ("gamma", "r--", False), ("Z", "g-.", False)]

fig, (ax1, ax2) = pyplot.subplots(1, 2)
analysis.plot_energy_distribution(ax1, Run, n_bins, resolution, Nsigma=Nsigma, comparisons=comparisons)
analysis.plot_energy_distribution(ax2, Run, n_bins, resolution, Nsigma=Nsigma, comparisons=[("tot", "k-", True)])
pyplot.suptitle("Monte Carlo scattering events: " + analysis.reaction_name)
if save_fig:
    size = fig.get_size_inches()
    size[0] = 1.75*size[0]
    size = size/scale
    fig.set_size_inches(size)
    pyplot.savefig(os.path.join(fig_path, f"energy_distribution.png"))
pyplot.show()

fig, ax = pyplot.subplots()
analysis.plot_energy_distribution(ax, Run, n_bins, resolution, theta_range=[0, numpy.pi/2], Nsigma=Nsigma, comparisons=comparisons)
pyplot.title("Forward Events: $\\theta \in [0, \\pi/2]$")
pyplot.show()

fig, ax = pyplot.subplots()
analysis.plot_energy_distribution(ax, Run, n_bins, resolution, theta_range=[numpy.pi/2, numpy.pi], Nsigma=Nsigma, comparisons=comparisons)
pyplot.title("Backward Events: $\\theta \in [\\pi/2, \\pi]$")
pyplot.show()

# ----------- FB asymmetry

fig, ax = pyplot.subplots()
analysis.plot_FB_asymmetry(ax, Run, n_bins, resolution, Nsigma)
if save_fig:
    size = fig.get_size_inches()
    size[0] = 1.75*size[0]
    size = size/scale
    fig.set_size_inches(size)
    pyplot.savefig(os.path.join(fig_path, f"FB_asymmetry.png"))
pyplot.show()
