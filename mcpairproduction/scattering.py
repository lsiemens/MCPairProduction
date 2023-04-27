"""Scattering cross sections and differential cross sections

The differential cross section and the analytic equation for the cross
section for components due to the photon the Z boson and the interference
of these two. These equations are derived for the two tree order diagrams
of the process e⁻ + e⁺ →  𝜇⁻ + 𝜇⁺. The analytic functions for the cross
section are given in terms of the indefinite integral of the differential
cross section so that the cross section over an arbitrary range of angles
can be found in addition to the total cross section. The differential
scattering cross sections are given in the center of momentum frame for
the ultrarelativistic limit and defined such that the angle θ is the
angle of the outgoing muon relative to the incoming electron.

While these equation where derived specifically for the process
e⁻ + e⁺ →  𝜇⁻ + 𝜇⁺ they are applicable more broadly. The equations can
be used for process of the form e⁻ + e⁺ →  f + f̄, where f is one of the
following fermions: muon, tauon, electron neutrino, muon neutrino and
tau neutrino. Note that in the case where the fermion is a neutrino then
only the Z boson component is relevant.

Note, internal calculations use natural units with ħ=c=1. Energies are
in units of MeV if not otherwise specified.

The differential cross sections for the photon and Z boson terms are
given in Griffiths [1]_ Section 8.1 Hadron production in e⁺e⁻ Collisions
and Example 9.5 Electron-Positron scattering near the Z pole. The full
equations for the differential cross section are given in section 2 of [2]_.

.. [1] Griffiths, David J. Introduction to Elementary Particles. 2nd,
       rev. ed. Weinheim: Wiley-VCH, 2008.

.. [2] Berends, F.A., R. Kleiss, and S. Jadach. "Radiative Corrections to
       Muon Pair and Quark Pair Production in Electron-Positron Collisions
       in the Z 0 Region." Nuclear physics. B 202, no. 1 (1982): 63-88.

"""

import numpy

from . import constants
from .constants import M_z, Gamma_z, g_e, g_z


def dsigma_dOmega(E, theta, A_2):
    """Differential cross section

    Compute the differential cross section in terms of the scattering
    amplitude squared. This equation is given for the center of momentum
    frame.

    The differential cross section dσ/dΩ is given by the equation below.

    dσ/dΩ = ⟨|A|²⟩/(16𝜋E)²

    Where the spin averaged scattering amplitude squared, ⟨|A|²⟩, is a
    function of the energy and the angle θ.

    Parameters
    ----------
    E : array or float
        The energy of the incoming particle in the center of momentum frame.
    theta :  array or float
        The angle, θ, between the incoming electron and the outgoing muon.
    A_2 : callable
        The spin averaged scattering amplitude squared.

    Returns
    -------
    dσ/dΩ : array or float
        The differential cross section of the process for the given
        energy and angle.
    """

    # the differential cross section
    return A_2(E, theta)/(16*numpy.pi*E)**2


def get_A_tot2(fermion_name="mu"):
    """Get ⟨|A|²⟩ function

    Get the total spin averaged scattering amplitude squared function
    for the interaction e⁻ + e⁺ →  f + f̄, where f is the specified
    fermion. The total amplitude squared can be expressed in terms of the
    following three components: photon term, Z boson term and the cross
    terms. The equation is given below,

    ⟨|A|²⟩ =  ⟨|A_𝛾|²⟩ + ⟨|A_Z|²⟩ + ⟨A_𝛾 A⁺_Z + A_Z A⁺_𝛾⟩

    Parameters
    ----------
    fermion_name : string
        Name of the resulting fundamental fermion in the reaction.
        Ex "mu" or "nu_e" ...

    Returns
    -------
    callable
        The function ⟨|A|²⟩ with arguments,
            - ``E``: The energy of the incoming particle.
            - ``theta``: The angle of the outgoing muon.
    """
    A_gamma2 = get_A_gamma2(fermion_name)
    A_Z2 = get_A_Z2(fermion_name)
    A_cross2 = get_A_cross2(fermion_name)

    def A_tot2(E, theta):
        """The function ⟨|A|²⟩

        The total spin averaged scattering amplitude squared.

        Parameters
        ----------
        E : array or float
            The energy of the incoming particle.
        theta : array or float
            The angle of the outgoing muon.
        """
        return A_gamma2(E, theta) + A_Z2(E, theta) + A_cross2(E, theta)
    return A_tot2


def get_A_gamma2(fermion_name="mu"):
    """Get ⟨|A_𝛾|²⟩ function

    Get the photon component of the spin averaged scattering amplitude
    squared function for the interaction e⁻ + e⁺ →  f + f̄, where f is
    the specified fermion. The equation is given below,

    ⟨|A_𝛾|²⟩ = g_e⁴*(1 + cos(θ)²)

    where g_e is the electromagnetic coupling constant.

    Parameters
    ----------
    fermion_name : string
        Name of the resulting fundamental fermion in the reaction.
        Ex "mu" or "nu_e" ...

    Returns
    -------
    callable
        The function ⟨|A_𝛾|²⟩ with arguments,
            - ``E``: The energy of the incoming particle.
            - ``theta``: The angle of the outgoing muon.
    """
    def A_gamma2(E, theta):
        """The function ⟨|A_𝛾|²⟩

        The photon component of the spin averaged scattering amplitude
        squared.

        Parameters
        ----------
        E : array or float
            The energy of the incoming particle.
        theta : array or float
            The angle of the outgoing muon.
        """
        return g_e**4*(1 + numpy.cos(theta)**2)
    return A_gamma2


def get_A_Z2(fermion_name="mu"):
    """Get ⟨|A_Z|²⟩ function

    Get the Z boson component of the spin averaged scattering amplitude
    squared function for the interaction e⁻ + e⁺ →  f + f̄, where f is
    the specified fermion. The equations are given below,

    C⁺ = (Cᵥ_e² + Cₐ_e²)(Cᵥ_𝜇² + Cₐ_𝜇²)
    Cˣ = Cᵥ_e Cₐ_e Cᵥ_𝜇 Cₐ_𝜇
    A = (g_z E)⁴/([4E² - M_z²]² + [M_z 𝚪_z]²)
    ⟨|A_Z|²⟩ = A (C⁺(1 + cos(θ)²) + 8Cˣcos(θ))

    where g_z is the neutral weak coupling, M_z is the mass of the Z boson,
    𝚪_z is the decay width of the Z boson and Cᵥ_e, Cₐ_e, Cᵥ_𝜇, Cₐ_𝜇 are
    the vector and axial-vector coupling constants of the electron and
    the other fermion.

    Parameters
    ----------
    fermion_name : string
        Name of the resulting fundamental fermion in the reaction.
        Ex "mu" or "nu_e" ...

    Returns
    -------
    callable
        The function ⟨|A_Z|²⟩ with arguments,
            - ``E``: The energy of the incoming particle.
            - ``theta``: The angle of the outgoing muon.
    """
    c_Ve, c_Ae, c_Vf, c_Af = get_neutral_couplings(fermion_name)
    c_sum = (c_Ve**2 + c_Ae**2)*(c_Vf**2 + c_Af**2)
    c_product = c_Ve*c_Ae*c_Vf*c_Af

    def A_Z2(E, theta):
        """The function ⟨|A_Z|²⟩

        The Z boson component of the spin averaged scattering amplitude
        squared.

        Parameters
        ----------
        E : array or float
            The energy of the incoming particle.
        theta : array or float
            The angle of the outgoing muon.
        """
        A = (g_z*E)**4/((4*E**2 - M_z**2)**2 + (M_z*Gamma_z)**2)
        return A*(c_sum*(1 + numpy.cos(theta)**2) + 8*c_product*numpy.cos(theta))
    return A_Z2


def get_A_cross2(fermion_name="mu"):
    """Get ⟨A_𝛾 A⁺_Z + A_Z A⁺_𝛾⟩ function

    Get the cross term component of the spin averaged scattering
    amplitude squared function for the interaction e⁻ + e⁺ →  f + f̄,
    where f is the specified fermion. The equations are given below,

    A = 8(g_e g_z E²)²/([4E² - M_z²]² + [M_z 𝚪_z]²)
    B = (1 - [M_z/(2 E)]²)
    ⟨A_𝛾 A⁺_Z + A_Z A⁺_𝛾⟩ = AB(Cᵥ_e Cᵥ_𝜇[1 + cos(θ)²] + 2Cₐ_e Cₐ_𝜇cos(θ))

    where g_e is the electromagnetic coupling, g_z is the neutral weak
    coupling, M_z is the mass of the Z boson, 𝚪_z is the decay width of
    the Z boson and Cᵥ_e, Cₐ_e, Cᵥ_𝜇, Cₐ_𝜇 are the vector and
    axial-vector coupling constants of the electron and the other fermion.

    Parameters
    ----------
    fermion_name : string
        Name of the resulting fundamental fermion in the reaction.
        Ex "mu" or "nu_e" ...

    Returns
    -------
    callable
        The function ⟨A_𝛾 A⁺_Z + A_Z A⁺_𝛾⟩ with arguments,
            - ``E``: The energy of the incoming particle.
            - ``theta``: The angle of the outgoing muon.
    """
    c_Ve, c_Ae, c_Vf, c_Af = get_neutral_couplings(fermion_name)

    def A_cross2(E, theta):
        """The function ⟨A_𝛾 A⁺_Z + A_Z A⁺_𝛾⟩

        The cross term component of the spin averaged scattering
        amplitude squared.

        Parameters
        ----------
        E : array or float
            The energy of the incoming particle.
        theta : array or float
            The angle of the outgoing muon.
        """
        A = 8*(g_e*g_z*E**2)**2/((4*E**2 - M_z**2)**2 + (M_z*Gamma_z)**2)
        B = (1 - (M_z/(2*E))**2)
        return A*B*(c_Ve*c_Vf*(1 + numpy.cos(theta)**2) + 2*c_Ae*c_Af*numpy.cos(theta))
    return A_cross2


def get_neutral_couplings(fermion_name):
    """Get coupling constants for the electron and other fermion

    Get the vector and axial-vector coupling constants for the electron
    and other specified fermion from `constants.py`. For more detail on
    the coupling constants check the docstring of
    `constants.neutral_couplings`.

    Parameters
    ----------
    fermion_name : string
        Name of the resulting fundamental fermion in the reaction.
        Ex "mu" or "nu_e" ...

    Returns
    -------
    (Cᵥ_e, Cₐ_e, Cᵥ_𝜇, Cₐ_𝜇) : tuple
        The vector and axial-vector couplings of the fermions.
    """
    if fermion_name == "e":
        raise ValueError(f"Invalid particle label: {fermion_name}\n"
                          "The reaction e⁻ + e⁺ →  e⁻ + e⁺ is not implemented.")

    # get the weak neutral coupling constants
    # coupling constants for the electron
    c_Ve, c_Ae = constants.neutral_couplings(*constants.fermions["e"])

    # fermion coupling constants for the other fermion
    c_Vf, c_Af = constants.neutral_couplings(*constants.fermions[fermion_name])

    return c_Ve, c_Ae, c_Vf, c_Af


def sigma_analytic(E, A_2_integrated):
    """Analytic cross section

    Compute the analytic cross section. The cross section is σ = ∬ (dσ/dΩ)dΩ
    where dσ/dΩ = ⟨|A|²⟩/(16𝜋E)². So expanding the integral then,
    σ = ∬ (dσ/dΩ)sin(θ)dθdφ = ∫⟨|A|²⟩sin(θ)dθ ∫(16𝜋E)⁻²dφ. Denoting the
    first term in this equation ⟨|A|²⟩_int then the analytic equation for
    the cross section is,

    σ = ⟨|A|²⟩_int/(128𝜋E²)

    Where ⟨|A|²⟩_int is a function of the energy.

    Parameters
    ----------
    E : array or float
        The energy of the incoming particle in the center of momentum frame.
    A_2_integrated : callable
        The integrated spin averaged scattering amplitude squared.

    Returns
    -------
    σ : array or float
        The cross section of the process for the given energy.
    """

    # the total cross section, integrated over phi
    return A_2_integrated(E)/(128*numpy.pi*E**2)


def get_A_tot2_integrated(fermion_name="mu", theta_range=[0, numpy.pi]):
    """Get ⟨|A|²⟩_int function

    Get the total integrated spin averaged scattering amplitude squared
    function for the interaction e⁻ + e⁺ →  f + f̄, where f is the
    specified fermion. The total amplitude squared can be expressed
    in terms of the following three components: photon term, Z boson term
    and the cross terms. The equation is given below,

    ⟨|A|²⟩_int = ⟨|A_𝛾|²⟩_int + ⟨|A_Z|²⟩_int + ⟨A_𝛾 A⁺_Z + A_Z A⁺_𝛾⟩_int

    where the terms are definite integrals of the form
    ⟨|A|²⟩_int = ∫⟨|A|²⟩sin(θ)dθ.

    Parameters
    ----------
    fermion_name : string
        Name of the resulting fundamental fermion in the reaction.
        Ex "mu" or "nu_e" ...
    theta_range : list
        The bounds of integration.

    Returns
    -------
    callable
        The function ⟨|A|²⟩_int with argument,
            - ``E``: The energy of the incoming particle.
    """
    A_gamma2_integrated = get_A_gamma2_integrated(fermion_name, theta_range)
    A_Z2_integrated = get_A_Z2_integrated(fermion_name, theta_range)
    A_cross2_integrated = get_A_cross2_integrated(fermion_name, theta_range)

    def A_tot2_integrated(E):
        """The function ⟨|A|²⟩_int

        The total integrated spin averaged scattering amplitude squared.

        Parameters
        ----------
        E : array or float
            The energy of the incoming particle.
        """
        return A_gamma2_integrated(E) + A_Z2_integrated(E) + A_cross2_integrated(E)
    return A_tot2_integrated


def get_A_gamma2_integrated(fermion_name="mu", theta_range=[0, numpy.pi]):
    """Get ⟨|A_𝛾|²⟩_int function

    Get the photon component of the integrated spin averaged scattering
    amplitude squared function for the interaction e⁻ + e⁺ →  f + f̄,
    where f is the specified fermion. The equation is given below,

    F(θ) = g_e⁴*(cos(θ) + cos(θ)³/3)
    ⟨|A_𝛾|²⟩_int = F(b) - F(a)

    where g_e is the electromagnetic coupling constant and [a, b] is the
    bounds of integration.

    Parameters
    ----------
    fermion_name : string
        Name of the resulting fundamental fermion in the reaction.
        Ex "mu" or "nu_e" ...
    theta_range : array or float
        The bounds of integration.

    Returns
    -------
    callable
        The function ⟨|A_𝛾|²⟩_int with argument,
            - ``E``: The energy of the incoming particle.
    """
    def A_gamma2_indefinite(E, theta):
        """The function ⟨|A_𝛾|²⟩_indefinite

        The indefinite integral of the photon component of the spin
        averaged scattering amplitude squared.

        Parameters
        ----------
        E : array or float
            The energy of the incoming particle.
        theta : array or float
            The angle of the outgoing muon.
        """
        return g_e**4*(numpy.cos(theta) + numpy.cos(theta)**3/3)
    return lambda E : A_gamma2_indefinite(E, theta_range[0]) - A_gamma2_indefinite(E, theta_range[1])


def get_A_Z2_integrated(fermion_name="mu", theta_range=[0, numpy.pi]):
    """Get ⟨|A_Z|²⟩_int function

    Get the Z boson component of the integrated spin averaged scattering
    amplitude squared function for the interaction e⁻ + e⁺ →  f + f̄,
    where f is the specified fermion. The equations are given below,

    C⁺ = (Cᵥ_e² + Cₐ_e²)(Cᵥ_𝜇² + Cₐ_𝜇²)
    Cˣ = Cᵥ_e Cₐ_e Cᵥ_𝜇 Cₐ_𝜇
    A = (g_z E)⁴/([4E² - M_z²]² + [M_z 𝚪_z]²)
    F(θ) = A (C⁺(cos(θ) + cos(θ)³/3) + 4Cˣcos(θ)²)
    ⟨|A_Z|²⟩ = F(b) - F(a)

    where g_z is the neutral weak coupling, M_z is the mass of the Z boson,
    𝚪_z is the decay width of the Z boson, [a, b] is the bounds of
    integration and Cᵥ_e, Cₐ_e, Cᵥ_𝜇, Cₐ_𝜇 are the vector and axial-vector
    coupling constants of the electron and the other fermion.

    Parameters
    ----------
    fermion_name : string
        Name of the resulting fundamental fermion in the reaction.
        Ex "mu" or "nu_e" ...
    theta_range : array or float
        The bounds of integration.

    Returns
    -------
    callable
        The function ⟨|A_Z|²⟩_int with argument,
            - ``E``: The energy of the incoming particle.
    """
    c_Ve, c_Ae, c_Vf, c_Af = get_neutral_couplings(fermion_name)
    c_sum = (c_Ve**2 + c_Ae**2)*(c_Vf**2 + c_Af**2)
    c_product = c_Ve*c_Ae*c_Vf*c_Af

    def A_Z2_indefinite(E, theta):
        """The function ⟨|A_Z|²⟩_indefinite

        The indefinite integral of the Z boson component of the spin
        averaged scattering amplitude squared.

        Parameters
        ----------
        E : array or float
            The energy of the incoming particle.
        theta : array or float
            The angle of the outgoing muon.
        """
        A = (g_z*E)**4/((4*E**2 - M_z**2)**2 + (M_z*Gamma_z)**2)
        return A*(c_sum*(numpy.cos(theta) + numpy.cos(theta)**3/3) + 4*c_product*numpy.cos(theta)**2)
    return lambda E : A_Z2_indefinite(E, theta_range[0]) - A_Z2_indefinite(E, theta_range[1])


def get_A_cross2_integrated(fermion_name="mu", theta_range=[0, numpy.pi]):
    """Get ⟨A_𝛾 A⁺_Z + A_Z A⁺_𝛾⟩_int function

    Get the cross term component of the integrated spin averaged scattering
    amplitude squared function for the interaction e⁻ + e⁺ →  f + f̄,
    where f is the specified fermion. The equations are given below,

    A = 8(g_e g_z E²)²/([4E² - M_z²]² + [M_z 𝚪_z]²)
    B = (1 - [M_z/(2 E)]²)
    F(θ) = AB(Cᵥ_e Cᵥ_𝜇[cos(θ) + cos(θ)³/3] + Cₐ_e Cₐ_𝜇cos(θ)²)
    ⟨A_𝛾 A⁺_Z + A_Z A⁺_𝛾⟩ = F(b) - F(a)

    where g_e is the electromagnetic coupling, g_z is the neutral weak
    coupling, M_z is the mass of the Z boson, 𝚪_z is the decay width of
    the Z boson, [a, b] is the bounds of integration and Cᵥ_e, Cₐ_e,
    Cᵥ_𝜇, Cₐ_𝜇 are the vector and axial-vector coupling constants of the
    electron and the other fermion.

    Parameters
    ----------
    fermion_name : string
        Name of the resulting fundamental fermion in the reaction.
        Ex "mu" or "nu_e" ...
    theta_range : array or float
        The bounds of integration.

    Returns
    -------
    callable
        The function ⟨A_𝛾 A⁺_Z + A_Z A⁺_𝛾⟩_int with argument,
            - ``E``: The energy of the incoming particle.
    """
    c_Ve, c_Ae, c_Vf, c_Af = get_neutral_couplings(fermion_name)

    def A_cross2_indefinite(E, theta):
        """The function  ⟨A_𝛾 A⁺_Z + A_Z A⁺_𝛾⟩_indefinite

        The indefinite integral of the cross term component of the spin
        averaged scattering amplitude squared.

        Parameters
        ----------
        E : array or float
            The energy of the incoming particle.
        theta : array or float
            The angle of the outgoing muon.
        """
        A = 8*(g_e*g_z*E**2)**2/((4*E**2 - M_z**2)**2 + (M_z*Gamma_z)**2)
        B = (1 - (M_z/(2*E))**2)
        return A*B*(c_Ve*c_Vf*(numpy.cos(theta) + numpy.cos(theta)**3/3) + c_Ae*c_Af*numpy.cos(theta)**2)
    return lambda E : A_cross2_indefinite(E, theta_range[0]) - A_cross2_indefinite(E, theta_range[1])
