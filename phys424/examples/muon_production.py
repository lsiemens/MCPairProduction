"""Estimat total cross section of muon production

Use totcross to estimate the total interaction cross section of the
neutral weak interaction e⁻ + e⁺ → 𝜇⁻ + 𝜇⁺.

"""

from matplotlib import pyplot
import numpy

from .. import totcross

def muon_Z_pole():
    # update latex preamble
    pyplot.rcParams.update({
        "font.family": "serif",
        "pgf.rcfonts": False,
        "pgf.texsystem": 'pdflatex'})

    M_z = totcross.M_z

    logM_z = numpy.log10(M_z/2)
    range_large = (3, 6)
    range_small = (85E3/2, 95E3/2)
    samples = 1000
    fermion_name = "mu"

    fig, (axleft, axright) = pyplot.subplots(1, 2)

    totcross.plot_compare(axleft, fermion_name, range_large, samples)
    totcross.plot_compare(axright, fermion_name, range_small, samples*20, logaxis=False)

    reaction_tex = totcross.get_reaction_equation(fermion_name)
    title = (f"Total cross section $\sigma(E)$ for {reaction_tex}"
              "\nMonte Carlo vs analytical integration")

    pyplot.suptitle(title)
    pyplot.show()
#    fig.set_size_inches(w=8, h=7)
#    pyplot.savefig("./scratch/Figure_1.pgf")
