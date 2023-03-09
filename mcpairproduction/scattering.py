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

from .constants import neutral_couplings, fermions, M_z, Gamma_z, g_z

def dsigma_dOmega(E, theta, fermion_name):
    """Differential cross section for electron positron scattering.

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
    # coupling constants for the electron
    c_Ve, c_Ae = neutral_couplings(*fermions["e"])

    # fermion coupling constants for the other fermion
    c_Vf, c_Af = neutral_couplings(*fermions[fermion_name])

    # precompute terms for the differential cross section
    A = (g_z**2*E/(16*numpy.pi))**2/((4*E**2 - M_z**2)**2 + (M_z*Gamma_z)**2)
    B = (c_Vf**2 + c_Af**2)*(c_Ve**2 + c_Ae**2)
    C = 8*c_Vf*c_Af*c_Ve*c_Ae

    # the differential cross section
    return A*(B*(1 + numpy.cos(theta)**2) - C*numpy.cos(theta))


def sigma_analytic(E, fermion_name):
    """Total cross section for electron positron scattering.

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
    # coupling constants for the electrons
    c_Ve, c_Ae = neutral_couplings(*fermions["e"])

    # fermion coupling constants for the other fermion
    c_Vf, c_Af = neutral_couplings(*fermions[fermion_name])

    # Calculate the total cross section
    A = (g_z**2*E/(16*numpy.pi))**2/((4*E**2 - M_z**2)**2 + (M_z*Gamma_z)**2)
    B = (c_Vf**2 + c_Af**2)*(c_Ve**2 + c_Ae**2)

    return 16*numpy.pi*A*B/3
