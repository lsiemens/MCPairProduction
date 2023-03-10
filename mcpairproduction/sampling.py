"""Sampling methods and Monte Carlo integration

The estimate of the total cross section is computed from the differential
cross section using Monte Carlo integration.

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

import numpy

def get_random_counts(expected_counts, seed=None):
    """Get random counts of events

    draw random event counts using poisson distribution
    """
    rng = numpy.random.default_rng(seed=seed)
    # currently (01/30/2023) default_rng uses the PCG-64 pseudo random
    # number generator by default.

    counts = rng.poisson(expected_counts)

    return counts

def get_random_samples(shape, seed=None):
    """Generate random samples of (θ, φ)

    The samples are generated uniformly in the domain, with θ ∊ [0, 𝜋)
    and φ ∊ [-𝜋, 𝜋). The size of the rectangular domain is also given.

    Parameters
    ----------
    shape : tuple
        Shape of array of sample point pairs. The resulting array of
        samples will have the shape (*shape, 2).
    seed : None, int
        Seed value for the random number generator. If `seed` is None,
        then numpy.random.default_rng will initalize with random entropy
        from the OS.

    Returns
    -------
    samples : array
        Random samples in the domain with the shape (*shape, 2).
    domain_size : float
        Size of the rectangular domain.
    """
    theta_min = 0
    theta_max = numpy.pi

    rng = numpy.random.default_rng(seed=seed)
    # currently (01/30/2023) default_rng uses the PCG-64 pseudo random
    # number generator by default.

    # TODO mention removing phi

    samples = rng.uniform(theta_min, theta_max, shape)
    # Note, random.Generator.uniform samples on a half open [low, high)
    # interval, so formaly samples with theta = 𝜋 or phi = 𝜋 are excluded.

    domain_size = 2*numpy.pi*numpy.prod(theta_max - theta_min)

    return samples, domain_size


def get_random_rejection_samples(shape, seed=None):
    """Generate random samples of (θ, v) for rejection sampling

    TODO mention removing phi

    The samples are generated uniformly in the domain, with θ ∊ [0, 𝜋),
    φ ∊ [-𝜋, 𝜋) and v ∊ [0, 1).

    Parameters
    ----------
    shape : tuple
        Shape of array of sample point pairs. The resulting array of
        samples will have the shape (2, *shape).
    seed : None, int
        Seed value for the random number generator. If `seed` is None,
        then numpy.random.default_rng will initalize with random entropy
        from the OS.

    Returns
    -------
    samples : array
        Random samples in the domain with the shape (2, *shape).
    """
    # remove phi [-pi, pi)
    theta_min_max = (0, numpy.pi)
    v_min_max = (0, 1)

    rng = numpy.random.default_rng(seed=seed)
    # currently (01/30/2023) default_rng uses the PCG-64 pseudo random
    # number generator by default.

    samples = numpy.empty(shape=(2, *shape))
    theta = rng.uniform(*theta_min_max, shape)
    v = rng.uniform(*v_min_max, shape)
    samples[0] = theta
    samples[1] = v
    # Note, random.Generator.uniform samples on a half open [low, high)
    # interval, so formaly samples with theta = 𝜋 or phi = 𝜋 are excluded.

    return samples


def monte_carlo_integration(function, x, dx, N):
    """MC integrate

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
    samples, domain_size = get_random_samples((N,))

    # Use monte carlo integration for each energy in E
    function_samples = function(x, samples)
    mean_sample = numpy.mean(function_samples*dx(samples))
    max_sample = numpy.max(function_samples)

    return mean_sample*domain_size, max_sample


def monte_carlo_integration_array(function, x, dx, N):
    """MC integrate

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
    samples, domain_size = get_random_samples((N, *x.shape))

    function_samples = function(x, samples)
    mean_sample = numpy.mean(function_samples*dx(samples), axis=0)
    max_sample = numpy.max(function_samples, axis=0)

    return mean_sample*domain_size, max_sample


def monte_carlo_sampling(function, x, maximum, N):
    """MC integrate

    Note x, maximum are not arrays
    """
    samples = numpy.empty(shape=(N,))
    missing = N
    missing_mask = numpy.ones(shape=(N,), dtype=numpy.bool)

    # recursively resample points that get rejected
    while missing > 0:
        new_samples_v = get_random_rejection_samples((missing, ))
        new_samples, new_v = new_samples_v[0], new_samples_v[1]
        new_values = function(x, new_samples)/maximum

        samples[missing_mask] = new_samples

        new_missing_mask = new_v > new_values
        new_missing = numpy.sum(new_missing_mask)
        missing_mask[missing_mask] = new_missing_mask
        missing = numpy.sum(missing_mask)
    return samples
