from matplotlib import pyplot
import numpy

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
