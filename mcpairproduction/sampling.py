"""Sampling methods and Monte Carlo integration

Rejection sampling and Monte Carlo integration are introduced in the PDG
book [1]_ section 42. For more on Monte Carlo integration and rejection
sampling is given by Gibbs [2]_ Chapter 2, see section 2.2 and 2.3.7

.. [1] R.L. Workman et al. (Particle Data Group), Prog. Theor. Exp. Phys.
   2022, 083C01 (2022)

.. [2] W.R. Gibbs, Computation in Modern Physics, 3rd Ed
"""

import numpy

def get_random_counts(expected_counts, seed=None):
    """Get random counts of events

    Draw random event counts using a Poisson distribution.

    Parameters
    ----------
    expected_counts : float
        The mean of the Poisson distribution.
    seed : None, integer
        Seed value for the random number generator. If `seed` is None,
        then numpy.random.default_rng will initialize with random entropy
        from the OS.

    Returns
    -------
    counts : integer
        A random integer sampled from the Poisson distribution.
    """
    rng = numpy.random.default_rng(seed=seed)
    # currently (01/30/2023) default_rng uses the PCG-64 pseudo random
    # number generator by default.

    counts = rng.poisson(expected_counts)

    return counts


def get_random_samples(shape, seed=None):
    """Generate random samples of an angle, θ

    The samples are generated uniformly in the domain, with θ ∊ [0, 𝜋).
    The size of the rectangular domain is also given.

    Parameters
    ----------
    shape : tuple
        Shape of array of sample points.
    seed : None, integer
        Seed value for the random number generator. If `seed` is None,
        then numpy.random.default_rng will initialize with random entropy
        from the OS.

    Returns
    -------
    samples : array
        Random samples in the domain.
    domain_size : float
        Size of the rectangular domain.
    """
    theta_min = 0
    theta_max = numpy.pi

    rng = numpy.random.default_rng(seed=seed)
    # currently (01/30/2023) default_rng uses the PCG-64 pseudo random
    # number generator by default.

    samples = rng.uniform(theta_min, theta_max, shape)
    # Note, random.Generator.uniform samples on a half open [low, high)
    # interval, so formaly samples with theta = 𝜋 are excluded.

    domain_size = 2*numpy.pi*(theta_max - theta_min)

    return samples, domain_size


def get_random_rejection_samples(shape, seed=None):
    """Generate random samples of (θ, v) for rejection sampling

    The samples are generated uniformly in the domain, with θ ∊ [0, 𝜋)
    and v ∊ [0, 1).

    Parameters
    ----------
    shape : tuple
        Shape of array of sample point pairs. The resulting array of
        samples will have the shape (2, *shape).
    seed : None, integer
        Seed value for the random number generator. If `seed` is None,
        then numpy.random.default_rng will initialize with random entropy
        from the OS.

    Returns
    -------
    samples : array
        Random samples in the domain with the shape (2, *shape).
    """
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
    # interval, so formaly samples with theta = 𝜋 or v = 1 are excluded.

    return samples


def dOmega(theta):
    """Solid angle element

    dΩ = sin(θ)dθdφ

    Parameters
    ----------
    theta : array or float
        Angular coordinate of the element

    """
    return numpy.sin(theta)


def monte_carlo_integration(function, fArg, N, dx=dOmega):
    """Monte Carlo integration

    An implementation of Monte Carlo integration on a spherical
    surface. Using `get_random_samples` to get N uniformly distributed
    random samples θ on the domain [0, 𝜋). Where the function, g(θ), to
    integrate over angle is function(θ)*dx(θ). The estimate of the
    integral is V*E[g(θ)], where E[] is the mean of the samples and V
    is the size of the domain. The maximum of g(θ) is computed for use
    in rejection sampling.

    Parameters
    ----------
    function : callable
        The function to integrate over the sphere.
    fArg : float
        Parameter to pass into `function`.
    N : integer
        The number of samples to use in the the Monte Carlo integral.
    dx : callable
        The differential element for the integral.

    Returns
    -------
    value : float
        Integral of `function` over all the sphere.
    max_sample : float
        The maximum sampled value of the integrand.
    """
    samples, domain_size = get_random_samples((N,))

    function_samples = function(fArg, samples)*dx(samples)
    mean_sample = numpy.mean(function_samples)
    max_sample = numpy.max(function_samples)

    return mean_sample*domain_size, max_sample


def monte_carlo_integration_array(function, fArg, N, dx=dOmega):
    """Monte Carlo integration for an array of parameters

    An implementation of Monte Carlo integration on a spherical
    surface. Using `get_random_samples` to get N uniformly distributed
    random samples θ on the domain [0, 𝜋). Where the function, g(θ), to
    integrate over angle is function(θ)*dx(θ). The estimate of the
    integral is V*E[g(θ)], where E[] is the mean of the samples and V
    is the size of the domain. The maximum of g(θ) is computed for use
    in rejection sampling.

    Parameters
    ----------
    function : callable
        The function to integrate over the sphere.
    fArg : array
        Parameter to pass into `function`.
    N : integer
        The number of samples to use in the the Monte Carlo integral.
    dx : callable
        The differential element for the integral.

    Returns
    -------
    value : array
        Integral of `function` over all the sphere.
    max_sample : array
        The maximum sampled value of the integrand.
    """
    samples, domain_size = get_random_samples((N, *fArg.shape))

    function_samples = function(fArg, samples)*dx(samples)
    mean_sample = numpy.mean(function_samples, axis=0)
    max_sample = numpy.max(function_samples, axis=0)

    return mean_sample*domain_size, max_sample


def monte_carlo_sampling(function, fArg, maximum, N, dx=dOmega):
    """Rejection sampling

    An implementation of rejection sampling on a spherical surface. Using
    `get_random_rejection_samples` to get N uniformly distributed random
    samples θ on the domain [0, 𝜋) along with corresponding random values
    v in the range [0, 1). Where the function, g(θ), to sample over angles
    is function(θ)*dx(θ). Any sample with v > g(θ)/g_max is rejected.
    This process is repeated to replace any rejected samples until all
    samples are accepted.

    Parameters
    ----------
    function : callable
        The function to sample on the sphere.
    fArg : float
        Parameter to pass into `function`.
    maximum : float
        The maximum of function(fArg, θ)dx(θ).
    N : integer
        The number of samples to generate.
    dx : callable
        The differential element

    Returns
    -------
    samples : array
        N random angles subject to the distribution `function(fArg, θ)`
        over the sphere.
    """
    samples = numpy.empty(shape=(N,))
    missing = N # number of samples that need to be generated
    # mask of samples that need to be generated
    missing_mask = numpy.ones(shape=(N,), dtype=numpy.bool)

    # recursively resample points that get rejected, untill no values
    # are missing (ie all samples are accepted).
    while missing > 0:
        # generate one sample for each one that is missing
        new_samples_v = get_random_rejection_samples((missing, ))
        new_samples, new_v = new_samples_v[0], new_samples_v[1]
        new_values = function(fArg, new_samples)*dx(new_samples)/maximum

        # use missing_mask to update samples with new_samples
        samples[missing_mask] = new_samples

        # mask of any of the new samples that are rejected
        new_missing_mask = new_v > new_values
        # count howmany on the new samples are rejected
        new_missing = numpy.sum(new_missing_mask)

        # use missing_mask to insert the new mask for the resampled values
        missing_mask[missing_mask] = new_missing_mask
        missing = numpy.sum(missing_mask)
    return samples
