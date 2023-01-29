"""
totcross.py: Estimate total cross sections

This script estimates the total interaction cross section for electron
positron scattering mediated by the Z boson e⁻ + e⁺ →  f + f̄. The
estimate of the total cross section is computed from the differential
cross section using Monte Carlo integration.

Note, internal calculations use natural units with ħ=c=1. Energies are
in units of MeV if not otherwise specified.

Refrence "Introduction to Elementary Particles" by Griffiths, 2nd Ed
Chapter 9 section 6: Neutral Weak Interactions

.. [1] D.J Griffiths "Introduction to Elementary Particles" 2nd Rev Ed. ...

.. [2] R.L. Workman et al. (Particle Data Group), Prog. Theor. Exp. Phys. 2022, 083C01 (2022)
"""

from matplotlib import pyplot
import numpy

# Constants: values from PDG
g_z = 0.7180 # TODO; the neutral coupling constant
sin2_theta_w = 0.2312 # unitless, sin²(θ_w) of the weak mixing angle θ_w
M_z = 9.118E4 # in MeV/c^2; mass of the Z boson
Gamma_z = 2.49E3 # in MeV/hbar; decay rate of Z boson

# a dictionary of fermion labels, all of the entries contain a tuple of
# the particle's weak isospin T_3 and charge Q. The tuples have the
# format (T_3, Q)
fermions = {"nu_e":( 0,  1/2), "nu_mu":( 0,  1/2), "nu_tau":( 0,  1/2),
               "e":(-1, -1/2),    "mu":(-1, -1/2),    "tau":(-1, -1/2),
            "u":( 2/3,  1/2), "c":( 2/3,  1/2), "t":( 2/3,  1/2),
            "d":(-1/3, -1/2), "s":(-1/3, -1/2), "b":(-1/3, -1/2)}

def neutral_couplings(T_3, Q):
    """Get the weak neutral couplings constants.

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
    """differential cross section for electron positron scattering

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

    TODO Add refrence to griffiths

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

#Note dΩ = sin(θ)dθdΦ

##### break: old code #####################################

#
# toy dsigmadOmega = (E/(E^2 - M_z^2))^2*(A(1 + cos^2(theta)) - Bcos(theta))
#
# degrees of freedom: theta, phi
#
# domin of variables: theta [0, pi]; phi [-pi, pi]
#

#def diff_cross(E, theta, phi):
#    return (E**2/((E**2 - M_z**2)**2 + Gamma**2))*(A*(1 + numpy.cos(theta)**2) - B*numpy.cos(theta))
#    return (A*(1 + numpy.cos(theta)**2) - B*numpy.cos(theta))

#M_z = 1
#A = 1
#B = 1
#Gamma = 0.01

N = 100

samples_theta = numpy.random.uniform(0, numpy.pi, size=N)
samples_phi = numpy.random.uniform(-numpy.pi, numpy.pi, size=N)

print((2*numpy.pi**2)*numpy.sum(numpy.sin(samples_theta)*dsigma_dOmega(M_z/2, samples_theta, "nu_e"))/N)

Es = numpy.linspace(0, M_z, 10*N)
thetas = numpy.linspace(0, numpy.pi, 10)
for theta in thetas:
    pyplot.plot(Es, dsigma_dOmega(Es, theta, "nu_e"))
pyplot.show()

Es = numpy.linspace(0, M_z, 10)
thetas = numpy.linspace(0, numpy.pi, 10*N)
for E in Es:
    pyplot.plot(thetas, dsigma_dOmega(E, thetas, "nu_e"))
pyplot.show()
