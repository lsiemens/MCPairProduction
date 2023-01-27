"""
totcross.py: Estimate total collision cross sections

Refrence "Introduction to Elementary Particles" by Griffiths, 2nd Ed
Chapter 9 section 6: Neutral Weak Interactions

In paticular compute the total cross section from the differential
cross section for e^+ + e^- -> f + \\bar{f}. Following Eample 9.5
"""

from matplotlib import pyplot
import numpy

# constants values from Griffiths
theta_w = 28.78 # in degrees, the weak mixing angle
sin2_theta_w = 0.2314 # unitless, sin^2(theta_w)

# todo select couplings from desired flaver
# Neutral vector and axial vector couplings
fermion_couplings = {"lepton":{""}, "quark":{}}


c_V = 1/2
c_A = 1/2
m_f = 0

#
# dsigma_dOmega = (g_z^2E/[16pi([2E]^2 - [M_z]^2)])^2 *
# {[(c_Vf)^2 + (c_Af)^2][(c_Ve)^2 + (c_Ae)^2](1 + cos^2(theta) - 8c_Vf c_Af c_Ve c_Ae cos(theta)}
#

##### break: old code #####################################

#
# toy dsigmadOmega = (E/(E^2 - M_z^2))^2*(A(1 + cos^2(theta)) - Bcos(theta))
#
# degrees of freedom: theta, phi
#
# domin of variables: theta [0, pi]; phi [-pi, pi]
#

def diff_cross(E, theta, phi):
#    return (E**2/((E**2 - M_z**2)**2 + Gamma**2))*(A*(1 + numpy.cos(theta)**2) - B*numpy.cos(theta))
    return (A*(1 + numpy.cos(theta)**2) - B*numpy.cos(theta))

M_z = 1
A = 1
B = 1
Gamma = 0.01

N = 100

samples_theta = numpy.random.uniform(0, numpy.pi, size=N)
samples_phi = numpy.random.uniform(-numpy.pi, numpy.pi, size=N)

print((2*numpy.pi**2)*numpy.sum(diff_cross(0.5, samples_theta, samples_phi))/N)

Es = numpy.linspace(0, 2, 1000)
thetas = numpy.linspace(0, numpy.pi, 10)
for theta in thetas:
    pyplot.plot(Es, diff_cross(Es, theta, 0))
pyplot.show()
