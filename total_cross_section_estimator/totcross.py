"""
totcross.py: Estimate total collision cross sections

Refrence "Introduction to Elementary Particles" by Griffiths, 2nd Ed
Chapter 9 section 6: Neutral Weak Interactions

In paticular compute the total cross section from the differential
cross section for e^+ + e^- -> f + \\bar{f}. Following Eample 9.5

Note internal calculations use natural units with h_bar=c=1.
Energies in units of MeV if not otherwise specified.
"""

from matplotlib import pyplot
import numpy

# constants values from Griffiths
g_z = 0.7180 # TODO; the neutral coupling constant
theta_w = 28.78 # in degrees, the weak mixing angle
sin2_theta_w = 0.2314 # unitless, sin^2(theta_w)
M_z = 9.1E4 # in MeV/c^2; mass of the Z boson
Gamma_z = 2.495E3 # in MeV/hbar; decay rate of Z boson

# todo select couplings from desired flaver
# Neutral vector and axial vector couplings
fermion_couplings = {"lepton":{""}, "quark":{}}

# neutrino
c_Vf = 1/2
c_Af = 1/2

# electron
c_Ve = -1/2 + 2*sin2_theta_w
c_Ae = -1/2

# TODO double check equation, in peticular the 16pi factor
# dsigma_dOmega = (g_z^2E)^2/([16pi[2E]^2 - [M_z]^2)]^2 + [16pi M_z Gamma_z]^2)
# *{ [(c_Vf)^2 + (c_Af)^2][(c_Ve)^2 + (c_Ae)^2](1 + cos^2(theta)
#    - 8c_Vf c_Af c_Ve c_Ae cos(theta)}
#

def dsigma_dOmega(E, theta):
    # breit-Wigner distribution
    BW_dist = (g_z**2*E/(16*numpy.pi))**2/((4*E**2 - M_z**2)**2 + (M_z*Gamma_z)**2)
    A = (c_Vf**2 + c_Af**2)*(c_Ve**2 + c_Vf**2)*(1 + numpy.cos(theta)**2)
    B = 8*c_Vf*c_Af*c_Ve*c_Ae*numpy.cos(theta)
    return BW_dist*(A - B)

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

print((2*numpy.pi**2)*numpy.sum(dsigma_dOmega(M_z/2, samples_theta))/N)

Es = numpy.linspace(0, M_z, 10*N)
thetas = numpy.linspace(0, numpy.pi, 10)
for theta in thetas:
    pyplot.plot(Es, dsigma_dOmega(Es, theta))
pyplot.show()

Es = numpy.linspace(0, M_z, 10)
thetas = numpy.linspace(0, numpy.pi, 10*N)
for E in Es:
    pyplot.plot(thetas, dsigma_dOmega(E, thetas))
pyplot.show()
