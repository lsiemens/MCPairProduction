# solotion for particle physics Assignment 1 question 2 part d
from matplotlib import pyplot
import numpy

latex_header = r"""
\documentclass{article}
\usepackage{graphicx}
\usepackage{fullpage}
\usepackage{natbib}
\usepackage{amsfonts}
\usepackage{amsmath}
\usepackage{hyperref}
\usepackage{xurl}
\usepackage{pgf}
\begin{document}"""

latex_footer = "\end{document}"

# particle formats: Name, Energy, Charge, Baryon, Strangeness
baryon_decuplet = [{"n":r"\Delta^-",	"e":1232, "c":-1, "s": 0},
                   {"n":r"\Delta^0",	"e":1232, "c": 0, "s": 0},
                   {"n":r"\Delta^+",	"e":1232, "c": 1, "s": 0},
                   {"n":r"\Delta^{++}", "e":1232, "c": 2, "s": 0},
                   {"n":r"\Sigma^{*-}",	"e":1387, "c":-1, "s":-1},
                   {"n":r"\Sigma^{*0}", "e":1384, "c": 0, "s":-1},
                   {"n":r"\Sigma^{*+}", "e":1382, "c": 1, "s":-1},
                   {"n":r"\Xi^{*-}",	"e":1535, "c":-1, "s":-2},
                   {"n":r"\Xi^{*0}",	"e":1535, "c": 0, "s":-2},
                   {"n":r"\Omega^{-}", 	"e":1672, "c":-1, "s":-3}]

baryon_octet = [{"n":r"n",	       "e": 940, "c": 0, "s": 0},
                {"n":r"p",         "e": 938, "c": 1, "s": 0},
                {"n":r"\Sigma^-",  "e":1197, "c":-1, "s":-1},
                {"n":r"\Sigma^0",  "e":1193, "c": 0, "s":-1},
                {"n":r"\Sigma^+",  "e":1189, "c": 1, "s":-1},
                {"n":r"\Lambda^0", "e":1116, "c": 0, "s":-1},
                {"n":r"\Xi^-",     "e":1322, "c":-1, "s":-2},
                {"n":r"\Xi^0",     "e":1315, "c": 0, "s":-2}]

meson_octet = [{"n":r"K^0",       "e": 498, "c": 0, "s": 1},
               {"n":r"K^+",       "e": 494, "c": 1, "s": 1},
               {"n":r"\pi^-",     "e": 140, "c":-1, "s": 0},
               {"n":r"\pi^0",     "e": 135, "c": 0, "s": 0},
               {"n":r"\pi^+",     "e": 140, "c": 1, "s": 0},
               {"n":r"\eta",      "e": 548, "c": 0, "s": 0},
               {"n":r"K^-",       "e": 594, "c":-1, "s":-1},
               {"n":"\\bar{K}^0", "e": 498, "c": 0, "s":-1}]

# List particles
body = "\\section{Baryons and Mesons}\n\\subsection{Baryon Decuplet}\n"
for baryon in baryon_decuplet:
    body += "$" + baryon["n"] + "$ "

body += "\n\n\\subsection{Baryon Octet}\n"

for baryon in baryon_octet:
    body += "$" + baryon["n"] + "$ "

body += "\n\n\\subsection{Meson Octet}\n"

for meson in meson_octet:
    body += "$" + meson["n"] + "$ "

# List reactions
body += "\n\n\\section{Reactions: quantum numbers only}\n"

for B_i in baryon_decuplet:
    for B_f in baryon_octet:
        for M_f in meson_octet:
            if (B_f["c"] + M_f["c"] - B_i["c"] == 0) and (B_f["s"] + M_f["s"] - B_i["s"] == 0):
                body += "$" + B_i["n"] + " \\rightarrow " + B_f["n"] + " + " + M_f["n"] + "$, "
    body += "\n\\newline\n"

text = latex_header + "\n\n" + body + "\n\n" + latex_footer
with open("./scratch/A1Q2d.tex", "w") as fout:
    fout.write(text)

