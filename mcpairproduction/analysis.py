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

from . import scattering
from . import constants
from . import utils

from matplotlib import pyplot
import numpy

reaction_name = f"$e^- + e^+ \\rightarrow \\mu^- + \\mu^+$"

differential_amplitudes = {"tot":(scattering.get_A_tot2, "$\\left|A\\right|^2$"),
                           "gamma":(scattering.get_A_gamma2, "$\\left|A_{\\gamma}\\right|^2$"),
                           "Z":(scattering.get_A_Z2, "$\\left|A_{Z}\\right|^2$"),
                           "cross":(scattering.get_A_cross2, "$A_{\\gamma}A_{Z}^* + A_{Z}A_{\\gamma}^*$")}

integrated_amplitudes = {"tot":(scattering.get_A_tot2_integrated, "$\\sigma_{total}$"),
                         "gamma":(scattering.get_A_gamma2_integrated, "$\\sigma_{\\gamma}$"),
                         "Z":(scattering.get_A_Z2_integrated, "$\\sigma_{Z}$"),
                         "cross":(scattering.get_A_cross2_integrated, "$\\sigma_{cross}$")}

def load(path):
    E, p_mag, theta, phi = numpy.loadtxt(path, delimiter=",", skiprows=2, unpack=True)
    with open(path, "r") as fin:
        next(fin)
        line = fin.readline()
        L_int = float(line)
    return E, p_mag, theta, phi, L_int


def setup_histogram(data, n_bins, range, Nsigma=1):
    hist, bin_edges = numpy.histogram(data, bins=n_bins, range=range)
    hist_error = Nsigma*numpy.sqrt(hist)
    bin_width = bin_edges[1] - bin_edges[0]

    return hist, hist_error, bin_edges[:-1], bin_width


def plot_angular_distribution(ax, Run_data, n_bins, resolution=100, Nsigma=1):
    E_data, p_data, theta_data, phi_data, L_int = Run_data

    E_data = E_data[0]
    theta_theory = numpy.linspace(0, numpy.pi, resolution)

    # calibrating the differential cross section to compair with the data
    tot_amplitude, _ = differential_amplitudes["tot"]
    ds_total = scattering.dsigma_dOmega(E_data, theta_theory, tot_amplitude())*numpy.sin(theta_theory)

    events_per_bin = len(theta_data)/n_bins
    normalization = events_per_bin/numpy.mean(ds_total)

    # make histogram of angular distribution
    hist, hist_err, bins, bin_width = setup_histogram(theta_data, n_bins, (0, numpy.pi), Nsigma)
    ax.bar(bins, hist, bin_width, align="edge", yerr=hist_err, color="b", ecolor="b", capsize=2, alpha=0.4, label="Data")

    sqrt_s = 2*E_data/1E3 # in GeV

    for component, style in [("tot", "k-"), ("gamma", "r--"), ("Z", "g-."), ("cross", "m:")]:
        diff_amplitude, label = differential_amplitudes[component]
        dsigma = scattering.dsigma_dOmega(E_data, theta_theory, diff_amplitude())*numpy.sin(theta_theory)
        ax.plot(theta_theory, normalization*dsigma, style, label=label)

    num_ticks = 4
    xticks = [i*numpy.pi/num_ticks for i in range(num_ticks + 1)]
    xtick_labels = [utils.nice_ticks(i, num_ticks) for i in range(num_ticks + 1)]

    ax.set_xticks(xticks, xtick_labels)
    ax.set_xlim(0, numpy.pi)
    ax.set_xlabel("$\\theta$ in radians")
    ax.set_ylabel("Scattering Events")
    ax.set_title(f"Angular distribution $\\frac{{d\\sigma}}{{d\\Omega}}$ at $\\sqrt{{s}} = {sqrt_s:.1f} GeV$")
    ax.legend()


def theoretic_binning(E, values, bins, bin_width):
    binned_values = []
    binned_E = []
    for bin in bins:
        mask = numpy.logical_and(E >= bin, E < bin + bin_width)
        binned_values.append(numpy.mean(values[mask]))
        binned_values.append(numpy.mean(values[mask]))
        binned_E.append(bin)
        binned_E.append(bin + bin_width)
    return binned_E, binned_values


def plot_energy_distribution(ax, Run_data, n_bins, resolution=100, theta_range=(0, numpy.pi), Nsigma=1, comparisons=[("tot", "k-", False)]):
    E_data, p_data, theta_data, phi_data, L_int = Run_data

    E_range = (numpy.min(E_data), numpy.max(E_data))
    E_theory = numpy.linspace(*E_range, resolution)

    mask = numpy.logical_and(theta_data >= theta_range[0], theta_data <= theta_range[1])
    hist, hist_err, bins, bin_width = setup_histogram(E_data[mask], n_bins, E_range, Nsigma)

    # convert to sqrt(s) in GeV
    bins, bin_width = 2*bins/1E3, 2*bin_width/1E3
    sqrt_s = 2*E_theory/1E3

#    ax.bar(bins, hist, bin_width, align="edge", yerr=hist_err, alpha=0.5, label="Data")
    ax.bar(bins, hist, bin_width, align="edge", yerr=hist_err, color="b", ecolor="b", capsize=2, alpha=0.4, label="Data")

    for component, style, binning in comparisons:
        int_amplitude, label = integrated_amplitudes[component]
        sigma = scattering.sigma_analytic(E_theory, int_amplitude(theta_range=theta_range))

        if binning:
            E_b, sigma_b = theoretic_binning(sqrt_s, sigma*L_int/n_bins, bins, bin_width)
            ax.semilogy(E_b, sigma_b, style, label=label + " binned")
        else:
            ax.semilogy(sqrt_s, sigma*L_int/n_bins, style, label=label)
    ax.set_ylim(1, None)
    ax.set_xlim(numpy.min(sqrt_s), numpy.max(sqrt_s))
    ax.set_ylabel("Scattering Events")
    ax.set_xlabel("$\\sqrt{s}$ in GeV")
    ax.legend()


def plot_FB_asymmetry(ax, Run_data, n_bins, resolution=100, Nsigma=1):
    E_data, p_data, theta_data, phi_data, L_int = Run_data

    E_range = (numpy.min(E_data), numpy.max(E_data))
    E_theory = numpy.linspace(*E_range, resolution)

    front_mask = theta_data <= numpy.pi/2
    back_mask = numpy.logical_not(front_mask)

    Fhist, Fhist_err, bins, bin_width = setup_histogram(E_data[front_mask], n_bins, E_range, Nsigma)
    Bhist, Bhist_err, bins, bin_width = setup_histogram(E_data[back_mask], n_bins, E_range, Nsigma)

    A_FB = (Fhist - Bhist)/(Fhist + Bhist)
    A_FB_err = (2*Fhist*Bhist/(Fhist + Bhist)**2)*numpy.sqrt((Fhist_err/Fhist)**2 + (Bhist_err/Bhist)**2)

    # convert to sqrt(s) in GeV
    bins, bin_width = 2*bins/1E3, 2*bin_width/1E3
    sqrt_s = 2*E_theory/1E3

#    ax.bar(bins, A_FB, bin_width, align="edge", yerr=A_FB_err, alpha=0.5, label="Data")
    ax.bar(bins, A_FB, bin_width, align="edge", yerr=A_FB_err, color="b", ecolor="b", capsize=2, alpha=0.4, label="Data")

    int_amplitude, _ = integrated_amplitudes["tot"]
    sigma_F = scattering.sigma_analytic(E_theory, int_amplitude(theta_range=[0, numpy.pi/2]))
    sigma_B = scattering.sigma_analytic(E_theory, int_amplitude(theta_range=[numpy.pi/2, numpy.pi]))
    FB_asymmetry = (sigma_F - sigma_B)/(sigma_F + sigma_B)

    ax.plot(sqrt_s, FB_asymmetry, "k-", label="$A_{FB}$")
    ax.set_xlim(numpy.min(sqrt_s), numpy.max(sqrt_s))
    ax.set_xlabel("$\\sqrt{s}$ in GeV")
    ax.set_ylabel("$A_{FB}$")
    ax.set_title("Forward Backward asymmetry: $A_{FB}$")
    ax.legend()

"""
def plot_compare(ax, fermion_name, range, MCsamples, resolution=100, logaxis=True):
    "
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
    "
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
"""
