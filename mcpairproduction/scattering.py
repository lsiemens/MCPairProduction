"""Scattering cross sections and differential cross sections

The differential cross section and the analytic equation for the cross
section for components due to the photon the Z boson and the interference
of these two. These equations are derived for the two tree order diagrams
of the process e⁻ + e⁺ →  𝜇⁻ + 𝜇⁺. The analytic functions for the cross
section are given interms of the indefinite integral of the differential
cross section so that the cross section over an arbitrary range of angles
can be found in adition to the total cross section. The differential
scattering cross sections are given in the centre of momentum frame for
the ultrarelitivistic limit and defined such that the angle θ is the
angle of the outgoing muon relative to the incoming electron.

While these equation where derived spesificaly for the process
e⁻ + e⁺ →  𝜇⁻ + 𝜇⁺ they are applicable more broadly. The equations can
be used for process of the form e⁻ + e⁺ →  f + f̄, where f is one of the
following fermions: muon, tauon, electron neutreno, muon neutreno and
tau neutreno. Note that in the case where the fermion is a neutrino then
only the Z boson componet is relevent.

Note, internal calculations use natural units with ħ=c=1. Energies are
in units of MeV if not otherwise specified.

The differentail cross sections for the photon and Z boson terms are
given in Griffiths [1]_ Section 8.1 Hadron production in e⁺e⁻ Collisions
and Example 9.5 Electron-Positron scattering near the Z pole. The full
equations for the differential cross section are given in section 2 of [2]_.

.. [1] Griffiths, David J. Introduction to Elementary Particles. 2nd,
       rev. ed. Weinheim: Wiley-VCH, 2008.

.. [2] Berends, F.A., R. Kleiss, and S. Jadach. “Radiative Corrections to
       Muon Pair and Quark Pair Production in Electron-Positron Collisions
       in the Z 0 Region.” Nuclear physics. B 202, no. 1 (1982): 63–88.

"""

import numpy

from . import constants
from .constants import M_z, Gamma_z, g_e, g_z

def dsigma_dOmega(E, theta, A_2):
    """Differential cross section

    Compute the differential cross section interms of the scattering
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
    """Get total ⟨|A|²⟩ function

    Get the total spin averaged scattering amplitude squared function
    for the interaction e⁻ + e⁺ →  f + f̄, where f is the specified
    fermion. The total amplitude squared can be expressed interms of the
    following three components: photon term, Z boson term and the
    interference term. The equation is given below,

    ⟨|A|²⟩ =  ⟨|A_𝛾|²⟩ + ⟨|A_Z|²⟩ + ⟨A_𝛾 A⁺_Z + A_Z A⁺_𝛾⟩

    Parameters
    ----------
    fermion_name : string
        Name of the resulting fundimental fermion in the reaction.
        Ex "mu" or "nu_e" ...

    Returns
    -------
    callable
        
        The differential cross section of the process for the given
        energy and angle.

    """
    A_gamma2 = get_A_gamma2(fermion_name)
    A_Z2 = get_A_Z2(fermion_name)
    A_cross2 = get_A_cross2(fermion_name)

    def A_tot2(E, theta):
        return A_gamma2(E, theta) + A_Z2(E, theta) + A_cross2(E, theta)
    return A_tot2

def get_A_gamma2(fermion_name="mu"):
    def A_gamma2(E, theta):
        return g_e**4*(1 + numpy.cos(theta)**2)
    return A_gamma2

def get_A_Z2(fermion_name="mu"):
    c_Ve, c_Ae, c_Vf, c_Af = get_neutral_couplings(fermion_name)
    c_sum = (c_Ve**2 + c_Ae**2)*(c_Vf**2 + c_Af**2)
    c_product = c_Ve*c_Ae*c_Vf*c_Af

    def A_Z2(E, theta):
        A = (g_z*E)**4/((4*E**2 - M_z**2)**2 + (M_z*Gamma_z)**2)
        return A*(c_sum*(1 + numpy.cos(theta)**2) + 8*c_product*numpy.cos(theta))
    return A_Z2

def get_A_cross2(fermion_name="mu"):
    """
    A = (g_z²E/[16𝜋])²/([4E² - M_z²]² + [M_z𝚪_z]²)
    B = (c_Vf² + c_Af²)(c_Ve² + c_Ae²)
    C = 8c_Vf c_Af c_Ve c_Ae

    """
    c_Ve, c_Ae, c_Vf, c_Af = get_neutral_couplings(fermion_name)

    def A_cross2(E, theta):
        A = 8*(g_e*g_z*E**2)**2/((4*E**2 - M_z**2)**2 + (M_z*Gamma_z)**2)
        B = (1 - (M_z/(2*E))**2)
        return A*B*(c_Ve*c_Vf*(1 + numpy.cos(theta)**2) + 2*c_Ae*c_Af*numpy.cos(theta))
    return A_cross2

def get_neutral_couplings(fermion_name):
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
    """Differential cross section in center of momentum frame and in
    the relitivistic limit.
    """

    # the total cross section, integrated over phi
    return A_2_integrated(E)/(128*numpy.pi*E**2)

def get_A_tot2_integrated(fermion_name="mu", theta_range=[0, numpy.pi]):
    A_gamma2_integrated = get_A_gamma2_integrated(fermion_name, theta_range)
    A_Z2_integrated = get_A_Z2_integrated(fermion_name, theta_range)
    A_cross2_integrated = get_A_cross2_integrated(fermion_name, theta_range)

    def A_tot2_integrated(E):
        return A_gamma2_integrated(E) + A_Z2_integrated(E) + A_cross2_integrated(E)
    return A_tot2_integrated

def get_A_gamma2_integrated(fermion_name="mu", theta_range=[0, numpy.pi]):
    def A_gamma2_indefinite(E, theta):
        return g_e**4*(numpy.cos(theta) + numpy.cos(theta)**3/3)
    return lambda E : A_gamma2_indefinite(E, theta_range[0]) - A_gamma2_indefinite(E, theta_range[1])

def get_A_Z2_integrated(fermion_name="mu", theta_range=[0, numpy.pi]):
    c_Ve, c_Ae, c_Vf, c_Af = get_neutral_couplings(fermion_name)
    c_sum = (c_Ve**2 + c_Ae**2)*(c_Vf**2 + c_Af**2)
    c_product = c_Ve*c_Ae*c_Vf*c_Af

    def A_Z2_indefinite(E, theta):
        A = (g_z*E)**4/((4*E**2 - M_z**2)**2 + (M_z*Gamma_z)**2)
        return A*(c_sum*(numpy.cos(theta) + numpy.cos(theta)**3/3) + 4*c_product*numpy.cos(theta)**2)
    return lambda E : A_Z2_indefinite(E, theta_range[0]) - A_Z2_indefinite(E, theta_range[1])

def get_A_cross2_integrated(fermion_name="mu", theta_range=[0, numpy.pi]):
    c_Ve, c_Ae, c_Vf, c_Af = get_neutral_couplings(fermion_name)

    def A_cross2_indefinite(E, theta):
        A = 8*(g_e*g_z*E**2)**2/((4*E**2 - M_z**2)**2 + (M_z*Gamma_z)**2)
        B = (1 - (M_z/(2*E))**2)
        return A*B*(c_Ve*c_Vf*(numpy.cos(theta) + numpy.cos(theta)**3/3) + c_Ae*c_Af*numpy.cos(theta)**2)
    return lambda E : A_cross2_indefinite(E, theta_range[0]) - A_cross2_indefinite(E, theta_range[1])
