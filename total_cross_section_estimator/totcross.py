"""
totcross.py: Estimate total cross sections

This script estimates the total interaction cross section for electron
positron scattering mediated by the Z boson e⁻ + e⁺ →  f + f̄. The
estimate of the total cross section is computed from the differential
cross section using Monte Carlo integration.

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

# Constants: values from PDG values
alpha = 7.29735256E-3 # unitless, fine-structure constant
sin2_theta_w = 0.2312 # unitless, sin²(θ_w) of the weak mixing angle θ_w
M_z = 9.118E4 # in MeV/c²; mass of the Z boson
Gamma_z = 2.49E3 # in MeV/hbar; decay rate of Z boson
hbarc2 = 3.893793721E5 # in MeV²mbarn

# calculating the coupling to the Z boson
# The electromagnetic coupling constant is g_e² = 4𝜋*𝛼 (note g_e equal
# to the elementary electric charge in natural units). This is related
# to the coupling constant for the Z boson g_z by the equation
# g_e = g_z cos(θ_w)sin(θ_w). Using these equations the constant g_z in
# terms of fine structure constant and weak mixing angle is
# g_z = sqrt(4𝜋*𝛼/[sin²(θ_w)(1 - sin²(θ_w))])
g_z = numpy.sqrt(4*numpy.pi*alpha/(sin2_theta_w - sin2_theta_w**2))

# a dictionary of fermion labels, all of the entries contain a tuple of
# the particle's weak isospin T_3 and charge Q. The tuples have the
# format (T_3, Q)
fermions = {"nu_e":( 0,  1/2), "nu_mu":( 0,  1/2), "nu_tau":( 0,  1/2),
               "e":(-1, -1/2),    "mu":(-1, -1/2),    "tau":(-1, -1/2),
            "u":( 2/3,  1/2), "c":( 2/3,  1/2), "t":( 2/3,  1/2),
            "d":(-1/3, -1/2), "s":(-1/3, -1/2), "b":(-1/3, -1/2)}

def neutral_couplings(T_3, Q):
    """
    Get the vector and axial-vector coupling constants.

    The coupling constants for a fermion are determined by the particle's
    weak isospin T₃ and charge Q, see PDG section 10.1 for more details.
    The vector coupling is,

    c_V = T₃ - 2Qsin²(θ_w)

    where θ_w is the weak mixing angle. The axial-vector coupling is,

    c_A = T₃

    Note in PDG material the coupling constants are reffered to as g_V and g_A.

    Parameters
    ----------
    T_3 : float
        The weak isospin of the fermion.
    Q : float
        The electric charge of the fermion.

    Returns
    -------
    (c_V, c_A) : tuple
        The vector and axial-vector coupling constants.
    """
    return T_3 - 2*Q*sin2_theta_w, T_3

def dsigma_dOmega(E, theta, fermion_name):
    """
    differential cross section for electron positron scattering

    Compute the differential cross section for the interaction
    e⁻ + e⁺ →  f + f̄ mediated by a Z boson. Note f may be any
    fundimental fermion other than the electron, that case has not
    been implemented.

    The differental cross section dσ/dΩ for these interactions is
    given by the equation below.

    A = (g_z²E/[16𝜋])²/([4E² - M_z²]² + [M_z𝚪_z]²)
    B = (c_Vf² + c_Af²)(c_Ve² + c_Ae²)
    C = 8c_Vf c_Af c_Ve c_Ae

    dσ/dΩ = A(B[1 + cos²(θ)] - C*cos(θ))

    where g_z is the neutral coupling constant, E is the energy of the
    incoming electron and positron, c_Ve, c_Ae are the neutral vector
    and axial vector couplings of the electron, c_Vf, c_Af are the
    neutral vector and axial vector couplings of resulting fermions,
    M_z is the mass of the Z boson and 𝚪_z is the decay rate of the Z boson.

    Parameters
    ----------
    E : array or float
        The energy of the incoming particles in the center of mass frame.
    theta : array or float
        The angle relative to the beam line of the resulting fermions.
    fermion_name : string
        Name of the resulting fundimental fermion for the reaction.
        Ex "mu" or "nu_e" ...

    Returns
    -------
    dσ/dΩ : array or float
        The differential cross section of the spesified interaction for
        the given E and theta.
    """

    if fermion_name == "e":
        raise ValueError(f"Invalid particle label: {fermion_name}\n"
                          "The reaction e⁻ + e⁺ →  e⁻ + e⁺ is not implemented.")

    # get the weak neutral coupling constants
    ## electron coupling constants
    c_Ve, c_Ae = neutral_couplings(*fermions["e"])

    ## fermion coupling constants
    c_Vf, c_Af = neutral_couplings(*fermions[fermion_name])

    # Calculate the differential cross section
    ## sub units
    A = (g_z**2*E/(16*numpy.pi))**2/((4*E**2 - M_z**2)**2 + (M_z*Gamma_z)**2)
    B = (c_Vf**2 + c_Af**2)*(c_Ve**2 + c_Vf**2)
    C = 8*c_Vf*c_Af*c_Ve*c_Ae

    ## the differential cross section
    return A*(B*(1 + numpy.cos(theta)**2) - C*numpy.cos(theta))

def sigma_analytic(E, fermion_name):
    """
    total cross section for electron positron scattering

    The total cross section computed by analytical integration of the
    differential cross section equation found in dsigma_dOmega.

    Parameters
    ----------
    E : array or float
        The energy of the incoming positron and electron.
    theta : array or float
        The angle relative to the beam line of the resulting fermions.
    fermion_name : string
        Name of the resulting fundimental fermion for the reaction.
        Ex "mu" or "nu_e" ...

    Returns
    -------
    σ : array or float
        The total cross section of the spesified interaction for the
        given E and theta.
    """
    # get the weak neutral coupling constants
    ## electron coupling constants
    c_Ve, c_Ae = neutral_couplings(*fermions["e"])

    ## fermion coupling constants
    c_Vf, c_Af = neutral_couplings(*fermions[fermion_name])

    # Calculate the total cross section
    A = (g_z**2*E/(16*numpy.pi))**2/((4*E**2 - M_z**2)**2 + (M_z*Gamma_z)**2)
    B = (c_Vf**2 + c_Af**2)*(c_Ve**2 + c_Vf**2)

    return 16*numpy.pi*A*B/3

def get_random_samples(shape, seed=None):
    """
    Generate random samples of (θ, φ)

    The samples are generated uniformly in the domain, with θ ∊ [0, 𝜋)
    and φ ∊ [-𝜋, 𝜋). The size of the rectangular domain is also given.

    Parameters
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
    theta_phi_min = numpy.array([0, -numpy.pi])
    theta_phi_max = numpy.array([numpy.pi, numpy.pi])

    rng = numpy.random.default_rng(seed=seed)
    # currently (01/30/2023) default_rng uses the PCG-64 pseudo random
    # number generator by default.

    samples = rng.uniform(theta_phi_min, theta_phi_max, (*shape, 2))
    # Note, random.Generator.uniform samples on a half open [low, high)
    # interval, so formaly samples with theta = 𝜋 or phi = 𝜋 are excluded.

    domain_size = numpy.prod(theta_phi_max - theta_phi_min)

    return samples, domain_size

def sigma_estimate(E, fermion_name, N):
    """
    Estimate the total cross section.

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
    """
    samples, domain_size = get_random_samples((*E.shape, N))

    theta = samples[:, :, 0]
    phi = samples[:, :, 1]

    # Use monte carlo integration for each energy in E
    function_samples = dsigma_dOmega(E[:, numpy.newaxis], theta, fermion_name)
    mean_sample = numpy.mean(function_samples*numpy.sin(theta), axis=1)
    #Note dΩ = sin(θ)dθdΦ

    return mean_sample*domain_size

def get_reaction_equation(fermion_name):
    """
    Get reaction equation for latex

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
    names = {"nu_e":("\\nu_e", "\\overline{\\nu}_e"),
             "nu_mu":("\\nu_\\mu", "\\overline{\\nu}_\\mu"),
             "nu_tau":("\\nu_\\tau", "\\overline{\\nu}_\\tau"),
             "e":("e^-", "e^+"),
             "mu":("\\mu^-", "\\mu^+"),
             "tau":("\\tau^-", "\\tau^+"),
             "u":("u", "\\overline{u}"), "d":("d", "\\overline{d}"),
             "c":("c", "\\overline{c}"), "s":("s", "\\overline{s}"),
             "t":("t", "\\overline{t}"), "b":("b", "\\overline{b}")}
    f, fbar = names[fermion_name]
    base_reaction = f"$e^+ + e^- \\rightarrow {f} + {fbar}$"
    return base_reaction

if __name__ == "__main__":
    # estimate the cross section on a logrithmic energy axis
    E_range_log = (4, 6)
    resolution = 100
    num_samples = 50

    # which reaction to use
    fermion_name = "mu"

    reaction_tex = get_reaction_equation(fermion_name)
    E_grid = 10**numpy.linspace(*E_range_log, resolution)
    E_grid_highres = 10**numpy.linspace(*E_range_log, 10*resolution)

    cross_section = hbarc2*sigma_estimate(E_grid, fermion_name, num_samples)
    analytic_cross_section = hbarc2*sigma_analytic(E_grid_highres, fermion_name)

    # graphing using 2*E_grid the total energy of the two particles
    pyplot.scatter(2*E_grid, cross_section, color="b", marker="+", label=f"$\\sigma_{{Monte Carlo}}$, {num_samples} samples")
    pyplot.loglog(2*E_grid_highres, analytic_cross_section, "k-", label="$\\sigma_{{Analytical}}$")

    pyplot.title(f"Total cross section $\sigma(E)$ for {reaction_tex}\nMonte Carlo vs analytical integration")
    pyplot.xlabel("Total Energy $E_{cm}$ in MeV")
    pyplot.ylabel("Total cross section in mbarn")
    pyplot.legend()
    pyplot.show()
