"""Validate new code

Reproduce the results of totcross using the new sampling and scattering
modules. Use the differential cross section for only the Z boson component
of the process e⁻ + e⁺ → 𝜇⁻ + 𝜇⁺ to validate the new code against the
the module totcross.py produced for the first part of the project.

"""

from matplotlib import pyplot
import numpy

from .. import totcross

from .. import sampling
from .. import scattering

fermion_name = "mu"

def dsigma_dOmega_Z(E, samples):
    """differential cross section using the new modules
    """
    A_Z2 = scattering.get_A_Z2(fermion_name)
    return scattering.dsigma_dOmega(E, samples, A_Z2)

def compare_energy_distributions():
    """Compare energy distrbitions

    Produce a plot to check if the new code can reproduce the energy
    distribution found by totcross. Note tot cross only uses the Z
    component of the differential cross section.
    """
    # Energy comparison muon
    M_z = totcross.M_z

    logM_z = numpy.log10(M_z/2)
    E_range = (85E3/2, 95E3/2)
    samples = 100000

    E = numpy.linspace(*E_range, 100)
    totcross_cross, _ = totcross.sigma_estimate(E, fermion_name, samples)
    sampling_cross, _ = sampling.monte_carlo_integration_array(dsigma_dOmega_Z, E, samples)
    analytic = totcross.sigma_analytic(E, fermion_name)
    pyplot.plot(E, totcross_cross, label="MC integral from totcross")
    pyplot.plot(E, sampling_cross, label="MC integral from sampling")
    pyplot.plot(E, analytic, label="MC analytic")
    pyplot.title("Comparison")
    pyplot.legend()
    pyplot.show()


def compare_angular_distributions():
    """Compare angular distrbitions

    Produce a plot to check if the new code can reproduce the angular
    distribution found by totcross. Note totcross only uses the Z
    component of the differential cross section.
    """
    # angular distribution comparison muon
    E = 45E3
    samples = 10000

    theta = numpy.linspace(0, numpy.pi, 1000)

    _, maximum = sampling.monte_carlo_integration(dsigma_dOmega_Z, E, samples)
    sampling_diffcross = sampling.monte_carlo_sampling(dsigma_dOmega_Z, E, maximum, samples)

    totcross_dist = totcross.dsigma_dOmega(E, theta, fermion_name)*numpy.sin(theta)
    pyplot.hist(sampling_diffcross, bins=50, alpha=0.5, density=True)
    pyplot.plot(theta, totcross_dist/numpy.mean(totcross_dist*numpy.pi))
    pyplot.show()
