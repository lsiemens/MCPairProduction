"""Physical constants and parameters

Internal calculations should use natural units with ħ=c=1. Energies are
in units of MeV if not otherwise specified. The PDG Book is the primary
source for the physical constants, *Review of Particle Physics* [1]_.

.. [1] R.L. Workman et al. (Particle Data Group), Prog. Theor. Exp. Phys.
   2022, 083C01 (2022)

"""

import numpy

# Leptons
M_e = 0.51099895000  # in MeV/c²; mass of the eletron
M_mu = 105.65837550  # in MeV/c²; mass of the muon
M_tau = 1776.86  # in MeV/c²; mass of the muon

# Bosons
M_w = 8.0377E4  # in MeV/c²; mass of the Z boson
Gamma_w = 2.085E3  # in MeV/hbar; decay rate of Z boson

M_z = 9.11876E4  # in MeV/c²; mass of the Z boson
Gamma_z = 2.4952E3  # in MeV/hbar; decay rate of Z boson

# Unit Conversions
hbarc2 = 3.893793721E5  # in MeV²mbarn

# unitless constants
alpha = 7.29735256E-3  # fine-structure constant
sin2_theta_w = 0.2312  # sin²(θ_w) where θ_w is the weak mixing angle

# calculating the coupling to the Z boson
# The electromagnetic coupling constant is g_e² = 4𝜋*𝛼 (note g_e equal
# to the elementary electric charge in natural units). This is related
# to the coupling constant for the Z boson g_z by the equation
# g_e = g_z cos(θ_w)sin(θ_w). Using these equations the constant g_z in
# terms of fine structure constant and weak mixing angle is
# g_z = sqrt(4𝜋*𝛼/[sin²(θ_w)(1 - sin²(θ_w))])
g_e = numpy.sqrt(4*numpy.pi*alpha)
g_z = numpy.sqrt(4*numpy.pi*alpha/(sin2_theta_w - sin2_theta_w**2))

# a dictionary of fermion labels, all of the entries contain a tuple of
# the particle's weak isospin T_3 and charge Q. The tuples have the
# format (T_3, Q)
fermions = {"nu_e": (1/2, 0),  "nu_mu": (1/2, 0), "nu_tau": (1/2, 0),
            "e": (-1/2, -1),   "mu": (-1/2, -1),  "tau": (-1/2, -1),
            "u": (1/2, 2/3),   "c": (1/2, 2/3),   "t": (1/2, 2/3),
            "d": (-1/2, -1/3), "s": (-1/2, -1/3), "b": (-1/2, -1/3)}


def neutral_couplings(T_3, Q):
    """Get the vector and axial-vector coupling constants

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
