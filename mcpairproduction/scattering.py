"""Scattering cross sections and differential cross sections

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

import numpy

from . import constants
from .constants import M_z, Gamma_z, g_e, g_z

def dsigma_dOmega(E, theta, A_2):
    """Differential cross section in center of momentum frame and in
    the relitivistic limit.
    """

    # the differential cross section
    return A_2(E, theta)/(16*numpy.pi*E)**2

def get_A_tot2(fermion_name="mu"):
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

def get_A_tot2_integrated(fermion_name="mu", theta_range=None):
    if theta_range is None:
        theta_range = [0, numpy.pi]

    A_gamma2_integrated = get_A_gamma2_integrated(fermion_name, theta_range)
    A_Z2_integrated = get_A_Z2_integrated(fermion_name, theta_range)
    A_cross2_integrated = get_A_cross2_integrated(fermion_name, theta_range)

    def A_tot2_integrated(E):
        return A_gamma2_integrated(E) + A_Z2_integrated(E) + A_cross2_integrated(E)
    return A_tot2_integrated

def get_A_gamma2_integrated(fermion_name="mu", theta_range=None):
    if theta_range is None:
        theta_range = [0, numpy.pi]

    def A_gamma2_indefinite(E, theta):
        return g_e**4*(numpy.cos(theta) + numpy.cos(theta)**3/3)
    return lambda E : A_gamma2_indefinite(E, theta_range[0]) - A_gamma2_indefinite(E, theta_range[1])

def get_A_Z2_integrated(fermion_name="mu", theta_range=None):
    if theta_range is None:
        theta_range = [0, numpy.pi]

    c_Ve, c_Ae, c_Vf, c_Af = get_neutral_couplings(fermion_name)
    c_sum = (c_Ve**2 + c_Ae**2)*(c_Vf**2 + c_Af**2)
    c_product = c_Ve*c_Ae*c_Vf*c_Af

    def A_Z2_indefinite(E, theta):
        A = (g_z*E)**4/((4*E**2 - M_z**2)**2 + (M_z*Gamma_z)**2)
        return A*(c_sum*(numpy.cos(theta) + numpy.cos(theta)**3/3) + 4*c_product*numpy.cos(theta)**2)
    return lambda E : A_Z2_indefinite(E, theta_range[0]) - A_Z2_indefinite(E, theta_range[1])

def get_A_cross2_integrated(fermion_name="mu", theta_range=None):
    if theta_range is None:
        theta_range = [0, numpy.pi]

    c_Ve, c_Ae, c_Vf, c_Af = get_neutral_couplings(fermion_name)

    def A_cross2_indefinite(E, theta):
        A = 8*(g_e*g_z*E**2)**2/((4*E**2 - M_z**2)**2 + (M_z*Gamma_z)**2)
        B = (1 - (M_z/(2*E))**2)
        return A*B*(c_Ve*c_Vf*(numpy.cos(theta) + numpy.cos(theta)**3/3) + c_Ae*c_Af*numpy.cos(theta)**2)
    return lambda E : A_cross2_indefinite(E, theta_range[0]) - A_cross2_indefinite(E, theta_range[1])
