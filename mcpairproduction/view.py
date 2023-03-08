"""view events

"""

import numpy
from matplotlib import pyplot

print("begin loading")
E, P, theta, phi = numpy.loadtxt("./scratch/out.dat", delimiter=",", unpack=True)
print("end loading")

n = theta.shape[0]
n_bins = 26
print(numpy.sqrt(n//n_bins))

pyplot.hist(theta, bins=n_bins)
pyplot.show()
pyplot.hist(phi, bins=n_bins)
pyplot.show()

hist, bin_edges = numpy.histogram(theta, bins=n_bins)

print(hist)
print(bin_edges)

A = hist
pyplot.bar(bin_edges[:-1], A, width=numpy.diff(bin_edges), yerr=numpy.sqrt(A))
pyplot.show()

A = (hist + hist[::-1])/2
pyplot.bar(bin_edges[:-1], A, width=numpy.diff(bin_edges), yerr=numpy.sqrt(A))
pyplot.show()

A = (hist - hist[::-1])/2
pyplot.bar(bin_edges[:-1], A, width=numpy.diff(bin_edges), yerr=numpy.sqrt(numpy.abs(A)))
pyplot.show()
