# solotion for particle physics Assignment 1 question 1 part f
from matplotlib import pyplot
import numpy

resolution = 1000

# assuming all masses   are in MeV/c^2
# assuming all momenta  are in MeV/c
# assuming all energies are in MeV

# mass of the electron
m_e = 0.5110 # in MeV/c^2

Q_3H    = 1.859E-2 # in MeV for    3H -> 3He + e-
Q_14C   = 1.565E-1 # in MeV for   14C -> 14N + e-
Q_137Cs = 1.176E+0 # in MeV for 137Cs -> 137Ba + e-

Qs = [Q_3H, Q_14C, Q_137Cs]
reactions = ["$\\ ^3H \\rightarrow \\ ^3He + e^+$",
             "$\\ ^{14}C \\rightarrow \\ ^{14}N + e^+$",
             "$\\ ^{137}Cs \\rightarrow \\ ^{137}Ba + e^+$"]

simple_name = ["3H_3He", "14C_14N", "137Cs_137Ba"]
line_styles = [("*k:", 80), ("sk-.", 100), ("ok--", 120)]


def find_zero(index, p, value):
    x1, y1 = p[index], value[index]
    x2, y2 = p[index + 1], value[index + 1]
    t = y1/(y1 - y2)
    interp = x1 + t*(x2 - x1)
    return interp

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

for Q, reaction, line_style in zip(Qs, reactions, line_styles):
    p = numpy.linspace(0, 1.1*p_m(Q), resolution)
    f_value = f(p,Q)
    pyplot.plot(p/p_m(Q), f_value/numpy.max(f_value), line_style[0], markevery=line_style[1], label=reaction)
pyplot.title("Comparison of f(p) on normalized axes")
pyplot.xlabel("relative momentum $p/p_c$ where $K(p_c) = Q$")
pyplot.ylabel("relative probability density $f(p) / f_{max}$")
pyplot.legend()
pyplot.savefig("./scratch/A1Q1_comp.png", dpi=300)
pyplot.show()

for Q, reaction, sname in zip(Qs, reactions, simple_name):
    p = numpy.linspace(0, 1.1*p_m(Q), resolution)
    f_value = f(p,Q)

    # normalize so the maximum is 1
    f_value = f_value/numpy.max(f_value)

    # find crossing points
    # subtract 0.5 so the value is zero when f(p) is at half the maximum
    f_value = f_value - 0.5

    # multiply the array by itself ofset by one element so the product
    # is negative if and only if the array crosses zero.
    cross = f_value[:-1]*f_value[1:]

    # set cross so the only non-zero values are where it is negative.
    # So cross is zero every where except wher the product of sequential
    # elements of f_value have oposing sign (ie crossing point of the
    # function) product.
    cross[cross > 0] = 0

    # find indices where cross is non-zero.
    indices = numpy.nonzero(cross)[0]

    # if "i" is an index where cross is negative, use linear interpolation
    # to find the point inbetween p[i] and p[i + 1] where the value of
    # f_value is zero (ie when the function is at half of it's maximum)
    i_min, i_max = indices
    p_min = find_zero(i_min, p, f_value)
    p_max = find_zero(i_max, p, f_value)
    p_range = numpy.array([p_min, p_max])
    p_half = numpy.array([p_range.mean(), p_max])

    HWHM = numpy.diff(p_range)[0]/2

    f_value = f(p,Q)
    f_max = numpy.max(f_value)
    print(f_max)

    hline = numpy.array([1/2, 1/2]) - 0.03

    pyplot.plot(p, f(p, Q)/f_max, "k-", label="$f(p)/f_{max}$")
    pyplot.plot(p_range, f(p_range, Q)/f_max, "k--")
    pyplot.plot(p_half, hline, "k:", label=f"HWHM: {HWHM:.2E} $MeV/c$")
    pyplot.scatter(p_half, hline, color="k", marker="|")
    pyplot.title("HWHM of $f(p)$ for " + reaction)
    pyplot.xlabel("momentum $p$ in $MeV/c$")
    pyplot.ylabel("relative probability density $f(p) / f_{max}$")
    pyplot.legend()
    pyplot.savefig(f"./scratch/A1Q1_{sname}.png", dpi=300)
    pyplot.show()
