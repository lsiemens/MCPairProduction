"""Analyze particle interaction event data

Estimate total cross sections of neutral weak interactions.

Estimates the total interaction cross section for electron positron
scattering mediated by the Z boson e⁻ + e⁺ →  f + f̄. The estimate of
the total cross section is computed from the differential cross section
using Monte Carlo integration.

Note, internal calculations use natural units with ħ=c=1. Energies are
in units of MeV if not otherwise specified.

The differentail cross section is from Griffiths [1]_ Example 9.5
Electron-Positron scattering near the Z pole, with *Review of Particle
Physics* [2]_ as an additional refrence and source for physical constants.
Equations for calculating the weak coupling to the Z boson are from
Thomson [3]_ section 15.3.1

.. [1] D. J. Griffiths, *Introduction to elementary particles*, 2nd,
   rev. ed. Weinheim: Wiley-VCH, 2008.

.. [2] R.L. Workman et al. (Particle Data Group), Prog. Theor. Exp. Phys.
   2022, 083C01 (2022)

.. [3] M. Thomson, *Modern particle physics*. Cambridge: Cambridge
   University Press, 2013.

"""

from matplotlib import pyplot
import numpy


def get_reaction_equation(fermion_name):
    """Get reaction equation for latex.

    Parameters
    ----------
    fermion_name : string
        Name of the resulting fundimental fermion for the reaction.
        Ex "mu" or "nu_e" ...

    Returns
    -------
    string
        Latex formated reaction equation.
    """
    names = {"nu_e": ("\\nu_e", "\\overline{\\nu}_e"),
             "nu_mu": ("\\nu_\\mu", "\\overline{\\nu}_\\mu"),
             "nu_tau": ("\\nu_\\tau", "\\overline{\\nu}_\\tau"),
             "e": ("e^-", "e^+"),
             "mu": ("\\mu^-", "\\mu^+"),
             "tau": ("\\tau^-", "\\tau^+"),
             "u": ("u", "\\overline{u}"), "d": ("d", "\\overline{d}"),
             "c": ("c", "\\overline{c}"), "s": ("s", "\\overline{s}"),
             "t": ("t", "\\overline{t}"), "b": ("b", "\\overline{b}")}
    f, fbar = names[fermion_name]
    base_reaction = f"$e^+ + e^- \\rightarrow {f} + {fbar}$"
    return base_reaction


def plot_compare(ax, fermion_name, range, MCsamples, resolution=100, logaxis=True):
    """
    Make plot of estimated total cross section

    Parameters
    ----------
    ax : pyplot axis
        The axis on which to make the plot.
    fermion_name : string
        Name of the resulting fundimental fermion for the reaction.
        Ex "mu" or "nu_e" ...
    range : tuple
        Range of energy to plot in MeV, or the logarithm of the energies
        if using logarithmic axes.
    MCsamples : int
        The number of samples to use in the the Monte Carlo integral for
        each energy value.
    resolution : int
        The resolution of sampled energy values.
    logaxis : boolean
        If true, the plot will be displaid with log-log axis.
    """
    if logaxis:
        E_MC = 10**numpy.linspace(*range, resolution)
        E_analytic = 10**numpy.linspace(*range, 10*resolution)
    else:
        E_MC = numpy.linspace(*range, resolution)
        E_analytic = numpy.linspace(*range, 10*resolution)

    cross_section, max_dsigma = sigma_estimate(E_MC, fermion_name, MCsamples)
    # convert cross_section from MeV² to mbarn
    cross_section = hbarc2*cross_section
    analytic = hbarc2*sigma_analytic(E_analytic, fermion_name)

    # graphing using the total energy of the two particles in GeV
    E_MC = 2*E_MC/1000
    E_analytic = 2*E_analytic/1000

    MC_label = f"$\\sigma_{{Monte\\ Carlo}}$, {MCsamples} samples"
    analytic_label = "$\\sigma_{{Analytical}}$"

    ax.scatter(E_MC, cross_section, color="b", marker="+", label=MC_label)
    if logaxis:
        ax.loglog(E_analytic, analytic, "k-", label=analytic_label)
    else:
        ax.plot(E_analytic, analytic, "k-", label=analytic_label)

    ax.plot(E_MC, hbarc2*4*numpy.pi*max_dsigma, label="max sample")

    ax.set_title(f"{E_MC[0]:.4g} GeV to {E_MC[-1]:.4g} GeV")
    ax.set_xlabel("Total Energy $E_{cm}$ in GeV")
    ax.set_ylabel("Total cross section in mbarn")
    ax.legend()
