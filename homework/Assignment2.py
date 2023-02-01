# Assignment 2
from matplotlib import pyplot
import numpy
import scipy.integrate

# Question 1.d
m_mu = 0.105659 # in GeV/c²
tau_rest = 2.2E-6 # in s, the rest frame mean life time of muons

def energy_dist(E, A=1):
    """
    Energy distribution of cosmic muons

    Using the energy distribution

    f(E) = AEe^-(E/σ)²

    with σ = 12/sqrt(𝜋) GeV

    Parameters
    ----------
    E : array
        Energy in GeV
    A : float
        Normalization constant in 1/GeV²
    """

    sigma = 12/numpy.sqrt(numpy.pi) # in GeV
    return A*E*numpy.exp(-(E/sigma)**2)

# the energy range from near the of the muon to the cutoff of 25 GeV
E_range = (m_mu + 1E-3, 25)
E = numpy.linspace(*E_range, 100)

A = 1/scipy.integrate.quad(energy_dist, *E_range)[0]

print(f"The normalization constant is {A} 1/GeV²")

expected_E = lambda E: E*energy_dist(E, A)

mean_E, _ = scipy.integrate.quad(expected_E, *E_range)
print(f"The mean energy is {mean_E:.6E} GeV")

pyplot.plot(E, energy_dist(E, A))
pyplot.show()

# Question 1.e
def beta(E):
    """
    relitivistic β factor for muons with energy E

    Parameters
    ----------
    E : array
        Energy in GeV
    """
    return numpy.sqrt(1 - (m_mu/E)**2)

def flux(f_0, d, E):
    """
    The flux some distance d towards the earth from the stratosphere
    d=0 km is an altitude of 50 km

    Parameters
    ----------
    f_0 : float
        The inital flux
    d : float
        distance traveled by the muons in km
    E : array
        Energy of the muons in GeV
    """

    c = 2.99792E5 # in km/s, the speed of light

    tau_earth = tau_rest/numpy.sqrt(1 - beta(E)**2)

    t = d/(c*beta(E))

    return f_0*numpy.exp(-t/tau_earth)

print(f"validate flux function\n\t{flux(1, -50, 6):.3f} should be ~3.88")
# diffrences due to truncation in the manual calculation

f_strat = 4 # in counts/(cm² minute)
expected_flux = lambda E: energy_dist(E, A)*flux(f_strat, 50, E)

mean_flux, _ = scipy.integrate.quad(expected_flux, *E_range)
print(f"The mean surface flux is {mean_flux} in counts/(cm² minute)")

expected_E = lambda E: E*expected_flux(E)/mean_flux
mean_E_surf, _ = scipy.integrate.quad(expected_E, *E_range)
print(f"The mean E of muons that arrive at the surface is {mean_E_surf} in GeV")
# as the muons pass deeper into the atmosphere the number that survice
# decreases but the average energy of those that servive increases. Not
# because any particles are gaining energy but the low energy ones are
# decaying first.

pyplot.plot(E, energy_dist(E, A)*f_strat)
pyplot.plot(E, expected_flux(E))
pyplot.show()
