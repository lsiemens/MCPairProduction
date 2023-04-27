#!/usr/bin/env python3
"""Analysis of LEP simulation data

This script produce plots of the simulated scattering events, and
optionally it prepares and saves the plots for the the slide presentation
on this project.
"""

from matplotlib import pyplot
import numpy
import os

from mcpairproduction import scattering
from mcpairproduction import analysis

if __name__ == "__main__":
    path = "./scratch/"

    # options for the slid presentation
    save_fig = False
    presentation = False
    if presentation:
        fig_path = "./presentation/Figures"
    else:
        fig_path = "../homework/phys424/project/Figures"
    scale = 1.4
    pyplot.rcParams.update({"savefig.bbox":"tight",
                            "savefig.dpi":250,
                            "font.family":"serif",
                            "pgf.rcfonts":False,
                            "pgf.texsystem":"pdflatex"})

    # Plot options
    resolution = 1000 # resolution of the theoretical curves
    n_bins = 32
    Nsigma = 2 # Size of error bars in sigma

    # read names of files produces during the last simulation run
    with open(os.path.join(path, "data_index"), "r") as fin:
        E_files = fin.read()
        E_files = E_files.split("\n")

    E_files, E_range_file = E_files[:-1], E_files[-1]

    Runs_fixed_E = [analysis.load(os.path.join(path, file)) for file in E_files]
    comparisons_angle = [("total", "k-", False), ("gamma", "r--", False),
                         ("Z", "g-.", False), ("cross", "m:", False)]

    # plots of the angular distribution at fixed energy
    for Run in Runs_fixed_E:
        fig, ax = pyplot.subplots()
        analysis.plot_angular_distribution(ax, Run, n_bins, resolution, Nsigma,
                                           comparisons=comparisons_angle)

        if save_fig:
            E, _, _, _, _ = Run
            E = E[0]
            fname = f"{E/1E3:.0f}GeV_angular_distribution"
            if presentation:
                size = fig.get_size_inches()
                size = size/scale
                fig.set_size_inches(size)
                ext = ".png"
                pyplot.savefig(os.path.join(fig_path, fname + ext))
            else:
                ext = ".pgf"
                pyplot.savefig(os.path.join(fig_path, fname + ext))
        pyplot.show()

    # plot the energy distribution for the energy sequence run and plots
    # of the forward portion of the events and likewise for the backward
    # portion of the events
    Run = analysis.load(os.path.join(path, E_range_file))
    comparisons_energy = [("total", "k-", False), ("gamma", "r--", False),
                          ("Z", "g-.", False)]

    fig, (ax1, ax2) = pyplot.subplots(1, 2)
    analysis.plot_energy_distribution(ax1, Run, n_bins, resolution, Nsigma=Nsigma,
                                      comparisons=comparisons_energy)
    analysis.plot_energy_distribution(ax2, Run, n_bins, resolution, Nsigma=Nsigma,
                                      comparisons=[("total", "k-", True)])

    pyplot.suptitle("Monte Carlo scattering events: " + analysis.reaction_name)
    if save_fig:
        if presentation:
            size = fig.get_size_inches()
            size[0] = 1.75*size[0]
            size = size/scale
            fig.set_size_inches(size)
            pyplot.savefig(os.path.join(fig_path, f"energy_distribution.png"))
        else:
            pyplot.savefig(os.path.join(fig_path, "energy_distribution.pgf"))
    pyplot.show()

    fig, ax = pyplot.subplots()
    analysis.plot_energy_distribution(ax, Run, n_bins, resolution,
                                      [0, numpy.pi/2], Nsigma,
                                      comparisons_energy)
    pyplot.title("Forward Events: $\\theta \in [0, \\pi/2]$")
    if save_fig and (not presentation):
        pyplot.savefig(os.path.join(fig_path, "forward.pgf"))
    pyplot.show()

    fig, ax = pyplot.subplots()
    analysis.plot_energy_distribution(ax, Run, n_bins, resolution,
                                      [numpy.pi/2, numpy.pi], Nsigma,
                                      comparisons_energy)
    pyplot.title("Backward Events: $\\theta \in [\\pi/2, \\pi]$")
    if save_fig and (not presentation):
        pyplot.savefig(os.path.join(fig_path, "backward.pgf"))
    pyplot.show()

    # plot the forward-backward asymmetry
    fig, ax = pyplot.subplots()
    analysis.plot_FB_asymmetry(ax, Run, n_bins, resolution, Nsigma)
    if save_fig:
        if presentation:
            size = fig.get_size_inches()
            size[0] = 1.75*size[0]
            size = size/scale
            fig.set_size_inches(size)
            pyplot.savefig(os.path.join(fig_path, "FB_asymmetry.png"))
        else:
            pyplot.savefig(os.path.join(fig_path, "FB_asymmetry.pgf"))
    pyplot.show()
