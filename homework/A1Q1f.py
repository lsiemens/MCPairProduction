# solotion for particle physics Assignment 1 question 1 part f
from matplotlib import pyplot
import numpy

resolution = 10

# assuming all masses   are in MeV/c^2
# assuming all momenta  are in MeV/c
# assuming all energies are in MeV

# mass of the electron
m_e = 0.511 # in MeV/c^2

Q_3H    = 1.859E-2 # in MeV for    3H -> 3He + e-
Q_14C   = 1.565E-1 # in MeV for   14C -> 14N + e-
Q_137Cs = 1.176E+0 # in MeV for 137Cs -> 137Ba + e-

Qs = [Q_3H, Q_14C, Q_137Cs]
reactions = ["$\\ ^3H \\rightarrow \\ ^3He + e^+$",
             "$\\ ^14C \\rightarrow \\ ^14N + e^+$",
             "$\\ ^137Cs \\rightarrow \\ ^137Ba + e^+$"]

def p_m(Q):
    p_max = numpy.sqrt(Q**2 + 2*m_e*Q)
    return p_max

def f(p, Q):

    E = numpy.sqrt(m_e**2 + p**2)
    f_value = p*(Q - E + m_e)**2*E

    # if p > p_m then set f_value to zero
    mask = (p > p_m(Q))
    f_value[mask] = 0

    return f_value

for Q, reaction in zip(Qs, reactions):
    p = numpy.linspace(0, 1.1*p_m(Q), resolution)
    print(Q, reaction)
    f_value = f(p,Q)
    pyplot.plot(p/p_m(Q), f_value/numpy.max(f_value), label=reaction)
pyplot.title("Comparison of f(p) on normalized axes")
pyplot.legend()
pyplot.show()

for Q, reaction in zip(Qs, reactions):
    p = numpy.linspace(0, 1.1*p_m(Q), resolution)
    print(Q, reaction)
    f_value = f(p,Q)
    # normalize so the maximum is 1
    f_value = f_value/numpy.max(f_value)

    # find crossing points
    # subtract 0.5 so the value is zero when f(p) is at half the maximum
    f_value = f_value - 0.5
    # multiply the array by itself ofset by one element so the product
    # is negative if and only if the array crosses zero.
    cross = f_value[:-1]*f_value[1:]
    # use the incises of these points to interpolate the crossing point
    # and the diffrence of the crossing points is the FWHM
    print(cross)
