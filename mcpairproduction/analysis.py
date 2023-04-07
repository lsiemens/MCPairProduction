"""Analyze particle scattering event data

This module contains functions to make plots of the scattering event
data. Plots of the angular distribution of events, plots of the energy
distribution of events and plots of the forward-backward asymmetry can
be produced. The forward-backward asymmetry characterizes the degree of 
asymmetry in the angular distribution at a given energy. It is defined as,

A_FB = (N_F - N_B)/(N_F + N_B)

where N_F is the number of events at a given energy with 0 < θ < 𝜋/2 and
N_B is the number of events at a given energy with 𝜋/2 < θ < 𝜋.

For more information about the forward-backward asymmetry see the PDG
book [1]_ section 10.5.1 Electro weak physics off the Z pole and Thomson
[3]_ chapter 16 section 2.2 Measurements of the weak mixing angle.

.. [1] R.L. Workman et al. (Particle Data Group), Prog. Theor. Exp. Phys.
   2022, 083C01 (2022)

.. [2] M. Thomson, *Modern particle physics*. Cambridge: Cambridge
   University Press, 2013.

"""

from . import scattering
from . import constants
from . import utils

from matplotlib import pyplot
import numpy

reaction_name = f"$e^- + e^+ \\rightarrow \\mu^- + \\mu^+$"

# Components of the spin averaged scattering amplitude squared
differential_amplitudes = {"total":(scattering.get_A_tot2, "$\\left|A\\right|^2$"),
                           "gamma":(scattering.get_A_gamma2, "$\\left|A_{\\gamma}\\right|^2$"),
                           "Z":(scattering.get_A_Z2, "$\\left|A_{Z}\\right|^2$"),
                           "cross":(scattering.get_A_cross2, "$A_{\\gamma}A_{Z}^* + A_{Z}A_{\\gamma}^*$")}

# Components of the integrated spin averaged scattering amplitude squared
integrated_amplitudes = {"total":(scattering.get_A_tot2_integrated, "$\\sigma_{total}$"),
                         "gamma":(scattering.get_A_gamma2_integrated, "$\\sigma_{\\gamma}$"),
                         "Z":(scattering.get_A_Z2_integrated, "$\\sigma_{Z}$"),
                         "cross":(scattering.get_A_cross2_integrated, "$\\sigma_{cross}$")}

def load(fname):
    """Load data

    Load scattering event data from produced by accelerator.py

    Parameters
    ----------
    fname : string
        Name of the data file.

    Returns
    -------
    E : array
        Energy of the particle.
    p_mag : array
        Magnitude of the 3-momentum of the particle.
    theta : array
        Angle of the 3-momentum relative to the incomming electron.
    phi : array
        Azimuthal angle of the particle, by convention the coordinate
        system is selected such that φ = 0.
    L_int : float
        The total integrated luminosity.
    """

    E, p_mag, theta, phi = numpy.loadtxt(fname, delimiter=",", skiprows=2, unpack=True)
    with open(fname, "r") as fin:
        next(fin)
        line = fin.readline()
        L_int = float(line)
    return E, p_mag, theta, phi, L_int


def setup_histogram(data, n_bins, range):
    """Create histogram

    Parameters
    ----------
    data : array
        Data to compute the histogram from.
    n_bins : integer
        The number of bins for the histogram of scattering events.
    range : tuple or list
        Range of the histogram.

    Returns
    -------
    hist : array
        Number of values in each bin.
    hist_err : array
        Error in the number of elements in each bin assuming Poisson
        statistics.
    bins : array
        Location of the left edge of each bin.
    bin_width : float
        Width of the histogram bins.
    """
    hist, bin_edges = numpy.histogram(data, bins=n_bins, range=range)
    hist_error = numpy.sqrt(hist)
    bin_width = bin_edges[1] - bin_edges[0]

    return hist, hist_error, bin_edges[:-1], bin_width


def theoretic_binning(x, y, bins, bin_width):
    """Bin theoretical curve

    Average the theoretic curve within each bin

    Note the resulting values for x, y produce a histogram like plot
    when plotted directly with pyplot.plot()

    Parameters
    ----------
    x : array
        x values of the theoretical curve.
    y : array
        y values of the theoretical curve.
    bins : array
        Location of the left edge of each bin.
    bin_width : float
        Width of the histogram bins.

    Returns
    -------
    binned_x : array
        Sequence of x values for each edge of each bin.
    binned_y
        Sequence of repeated y values. The mean value of y in the bin,
        for each edge of each bin.
    """
    binned_x = []
    binned_y = []
    for bin in bins:
        mask = numpy.logical_and(x >= bin, x < bin + bin_width)
        mean_y = numpy.mean(y[mask])
        binned_x = binned_x + [bin, bin + bin_width]
        binned_y = binned_y + [mean_y, mean_y]
    return binned_x, binned_y


def plot_angular_distribution(ax, Run_data, n_bins, resolution=100, Nsigma=1, comparisons=[("total", "k-", False)]):
    """Plot scattering events vs angle

    Make plot of scattering events vs angle and compare against theoretical
    predictions of the expected number of events.

    Parameters
    ----------
    ax : pyplot axis
        The axis on which to make the plot.
    Run_data : tuple
        Tuple of data from loading a data set using `analysis.load`.
    n_bins : integer
        The number of bins for the histogram of scattering events.
    resolution : int
        The resolution for theoretical curves.
    Nsigma : integer
        show error bars with size `Nsigma` standard deviations.
    comparisons : list
        Set of theoretical curves for comparison where each curve is given
        as (Name, line_style, bin_curve). Name is one of
        ["total", "gamma", "Z", "cross"], line_style is a matplotlib line
        style string and if bin_curve is True then the theoretical curve
        is averaged over each bin.
    """
    E_data, p_data, theta_data, phi_data, L_int = Run_data

    E_data = E_data[0]
    theta_theory = numpy.linspace(0, numpy.pi, resolution)

    # calibrating the differential cross section to compair with the data
    tot_amplitude, _ = differential_amplitudes["total"]
    ds_total = scattering.dsigma_dOmega(E_data, theta_theory, tot_amplitude())*numpy.sin(theta_theory)

    events_per_bin = len(theta_data)/n_bins
    normalization = events_per_bin/numpy.mean(ds_total)

    # make histogram of angular distribution
    hist, hist_err, bins, bin_width = setup_histogram(theta_data, n_bins, (0, numpy.pi))
    ax.bar(bins, hist, bin_width, align="edge", yerr=Nsigma*hist_err, color="b", ecolor="b", capsize=2, alpha=0.4, label="Data")

    sqrt_s = 2*E_data/1E3 # in GeV

    for component, style, binning in comparisons:
        diff_amplitude, label = differential_amplitudes[component]
        dsigma = scattering.dsigma_dOmega(E_data, theta_theory, diff_amplitude())*numpy.sin(theta_theory)

        if binning:
            theta_b, dsigma_b = theoretic_binning(theta_theory, normalization*dsigma, bins, bin_width)
            ax.plot(theta_b, dsigma_b, style, label=label)
        else:
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


def plot_energy_distribution(ax, Run_data, n_bins, resolution=100, theta_range=(0, numpy.pi), Nsigma=1, comparisons=[("total", "k-", False)]):
    """Plot scattering events vs energy

    Make plot of scattering events vs energy and compare against theoretical
    predictions of the expected number of events.

    Parameters
    ----------
    ax : pyplot axis
        The axis on which to make the plot.
    Run_data : tuple
        Tuple of data from loading a data set using `analysis.load`.
    n_bins : integer
        The number of bins for the histogram of scattering events.
    resolution : int
        The resolution for theoretical curves.
    theta_range : tuple
        Consider only events within this range of angles.
    Nsigma : integer
        show error bars with size `Nsigma` standard deviations.
    comparisons : list
        Set of theoretical curves for comparison where each curve is given
        as (Name, line_style, bin_curve). Name is one of
        ["total", "gamma", "Z", "cross"], line_style is a matplotlib line
        style string and if bin_curve is True then the theoretical curve
        is averaged over each bin.
    """

    E_data, p_data, theta_data, phi_data, L_int = Run_data

    E_range = (numpy.min(E_data), numpy.max(E_data))
    E_theory = numpy.linspace(*E_range, resolution)

    mask = numpy.logical_and(theta_data >= theta_range[0], theta_data <= theta_range[1])
    hist, hist_err, bins, bin_width = setup_histogram(E_data[mask], n_bins, E_range)

    # convert to sqrt(s) in GeV
    bins, bin_width = 2*bins/1E3, 2*bin_width/1E3
    sqrt_s = 2*E_theory/1E3

    ax.bar(bins, hist, bin_width, align="edge", yerr=Nsigma*hist_err, color="b", ecolor="b", capsize=2, alpha=0.4, label="Data")

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
    """Plot the Forward-Backward asymmetry

    Make plot of forward-backward asymetry vs energy and compare against
    theoretical predictions.

    Parameters
    ----------
    ax : pyplot axis
        The axis on which to make the plot.
    Run_data : tuple
        Tuple of data from loading a data set using `analysis.load`.
    n_bins : integer
        The number of bins for the histogram of scattering events.
    resolution : int
        The resolution for theoretical curves.
    Nsigma : integer
        show error bars with size `Nsigma` standard deviations.
    """
    E_data, p_data, theta_data, phi_data, L_int = Run_data

    E_range = (numpy.min(E_data), numpy.max(E_data))
    E_theory = numpy.linspace(*E_range, resolution)

    front_mask = theta_data <= numpy.pi/2
    back_mask = numpy.logical_not(front_mask)

    Fhist, Fhist_err, bins, bin_width = setup_histogram(E_data[front_mask], n_bins, E_range)
    Bhist, Bhist_err, bins, bin_width = setup_histogram(E_data[back_mask], n_bins, E_range)

    A_FB = (Fhist - Bhist)/(Fhist + Bhist)
    A_FB_err = (2*Fhist*Bhist/(Fhist + Bhist)**2)*numpy.sqrt((Fhist_err/Fhist)**2 + (Bhist_err/Bhist)**2)

    # convert to sqrt(s) in GeV
    bins, bin_width = 2*bins/1E3, 2*bin_width/1E3
    sqrt_s = 2*E_theory/1E3

    ax.bar(bins, A_FB, bin_width, align="edge", yerr=Nsigma*A_FB_err, color="b", ecolor="b", capsize=2, alpha=0.4, label="Data")

    int_amplitude, _ = integrated_amplitudes["total"]
    sigma_F = scattering.sigma_analytic(E_theory, int_amplitude(theta_range=[0, numpy.pi/2]))
    sigma_B = scattering.sigma_analytic(E_theory, int_amplitude(theta_range=[numpy.pi/2, numpy.pi]))
    FB_asymmetry = (sigma_F - sigma_B)/(sigma_F + sigma_B)

    ax.plot(sqrt_s, FB_asymmetry, "k-", label="$A_{FB}$")
    ax.set_xlim(numpy.min(sqrt_s), numpy.max(sqrt_s))
    ax.set_xlabel("$\\sqrt{s}$ in GeV")
    ax.set_ylabel("$A_{FB}$")
    ax.set_title("Forward Backward asymmetry: $A_{FB}$")
    ax.legend()
