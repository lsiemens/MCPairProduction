""" luminosity / write data / ...

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

def sigma_estimate(E, fermion_name, N):
    """Estimate the total cross section.

    Use Monte Carlo integration to find the total cross section from the
    differential cross section dσ/dΩ, note dΩ = sin(θ)dθdφ. Then the
    total cross section is σ = ∬ (dσ/dΩ)sin(θ)dθdφ over the rectangular
    domain with θ ∊ [0, 𝜋) and φ ∊ [-𝜋, 𝜋).

    Parameters
    ----------
    E : array
        An array of energies to evaluat the total cross section at.
    fermion_name : string
        Name of the resulting fundimental fermion for the reaction.
        Ex "mu" or "nu_e" ...
    N : int
        The number of samples to use in the the Monte Carlo integral for
        each energy value.

    Returns
    -------
    σ : array
        An array of estimates of the total cross section with the same
        shape as the array `E`.
    max_sample : float
        The maximum sampled value of dσ/dΩ
    """
    samples, domain_size = get_random_samples((*E.shape, N))

    theta = samples[:, :, 0]
    phi = samples[:, :, 1]

    # Use monte carlo integration for each energy in E
    function_samples = dsigma_dOmega(E[:, numpy.newaxis], theta, fermion_name)
    mean_sample = numpy.mean(function_samples*numpy.sin(theta), axis=1)
    max_sample = numpy.max(function_samples, axis=1)
    # Note dΩ = sin(θ)dθdφ

    return mean_sample*domain_size, max_sample
