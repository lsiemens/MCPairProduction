"""tau mc event

"""

import numpy
from matplotlib import pyplot
from . import totcross

E = 90E3 # MeV
E = numpy.array([E])
N = 10000
m_tau = 1776.86 # MeV

t = numpy.linspace(0, numpy.pi, 100)

totsigma, max = totcross.sigma_estimate(E, "tau", 1000)
print(totsigma, max, "max")

samples = totcross.get_random_rejection_samples((N,))
theta = samples[:, 0]
phi = samples[:, 1]
x = samples[:, 2]

values = totcross.dsigma_dOmega(E, theta, "tau")/max

map = values >= x
nmap = numpy.logical_not(map)

pyplot.plot(t, totcross.dsigma_dOmega(E, t, "tau")/max, "k")
pyplot.scatter(theta[map], x[map])
pyplot.scatter(theta[nmap], x[nmap])
pyplot.show()

print("Begin saving")
with open("./scratch/out.dat", "w") as fout:
    string = "# Data: E, P, theta, phi\n"
    string += "\n".join([f"{E[0]:12.5E},{numpy.sqrt(E[0]**2 - m_tau**2):12.5E},{th:12.5E},{ph:12.5E}" for th, ph in zip(theta[map], phi[map])])
    fout.write(string)
print("end saving")

