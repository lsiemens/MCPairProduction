"""Generate letex documents with Feynman digrams using feynMP
"""

import os
import subprocess

class Environment:
    def __init__(self, type, optional=None):
        self.type=type
        self.contents = []
        self.optional = optional

    def __str__(self):
        string = f"\\begin{{{self.type}}}"
        if self.optional is not None:
            string += str(self.optional)
        string += "\n"

        string += "\n".join([str(content) for content in self.contents])
        string += f"\n\\end{{{self.type}}}"
        return string

    def add(self, content):
        self.contents.append(content)

class Document:
    def __init__(self):
        self.dclass = "article"
        self.packages = ["feynmp"]

        header = "\\DeclareGraphicsRule{*}{mps}{*}{}"
        self.contents = [header]

    def write(self, path):
        with open(path, "w") as fin:
            fin.write(str(self))

    def compile(self, path, mp_name, display=False):
        cwd, file = os.path.split(path)
        subprocess.run(["pdflatex", file], cwd=cwd)
        subprocess.run(["mpost", mp_name], cwd=cwd)
        subprocess.run(["pdflatex", file], cwd=cwd)

        if display:
            pdf_file, _ = os.path.splitext(path)
            subprocess.run(["xdg-open", pdf_file + ".pdf"])

    def add(self, content):
        self.contents.append(content)

    def __str__(self):
        string = f"\\documentclass{{{self.dclass}}}\n"
        string += "\n".join([f"\\usepackage{{{package}}}" for package in self.packages])
        string += "\n"

        string += "\n".join([str(content) for content in self.contents])
        return string

def left(vertices):
    string = f"\\fmfleft{{{','.join(vertices)}}}"
    return string

def right(vertices):
    string = f"\\fmfright{{{','.join(vertices)}}}"
    return string

def fermion(vertices):
    string = f"\\fmf{{fermion}}{{{','.join(vertices)}}}"
    return string

def photon(vertices):
    string = f"\\fmf{{photon}}{{{','.join(vertices)}}}"
    return string

test = Document()

doc = Environment("document")
test.add(doc)
test_mp = Environment("fmffile", "{test_mp}")
doc.add(test_mp)
graph = Environment("fmfgraph", "(120, 80)")
test_mp.add(graph)
graph.add(left(["i1", "i2"]))
graph.add(right(["o1", "o2"]))
graph.add(fermion(["i1", "v1", "i2"]))
graph.add(fermion(["o1", "v2", "o2"]))
graph.add(photon(["v1", "v2"]))

graph2 = Environment("fmfgraph", "(120, 80)")
test_mp.add(graph2)
graph2.add(left(["i1", "i2"]))
graph2.add(right(["o1", "o2"]))
graph2.add(fermion(["i1", "v1", "o1"]))
graph2.add(fermion(["i2", "v2", "o2"]))
graph2.add(photon(["v1", "v2"]))

graph3 = Environment("fmfgraph", "(120, 80)")
test_mp.add(graph3)
graph3.add(left(["i1", "i2"]))
graph3.add(right(["o1", "o2"]))
graph3.add(fermion(["i1", "v1", "v2", "o1"]))
graph3.add(fermion(["i2", "v3", "v4", "o2"]))
graph3.add(photon(["v1", "v3"]))
graph3.add(photon(["v2", "v4"]))

graph4 = Environment("fmfgraph", "(120, 80)")
test_mp.add(graph4)
graph4.add(left(["i1", "i2"]))
graph4.add(right(["o1", "o2"]))
graph4.add(fermion(["i1", "v1", "v2", "i2"]))
graph4.add(fermion(["o1", "v3", "v4", "o2"]))
graph4.add(photon(["v1", "v3"]))
graph4.add(photon(["v2", "v4"]))

# problems with swaping external particles
graph6 = Environment("fmfgraph", "(120, 80)")
test_mp.add(graph6)
graph6.add(left(["i1", "i2"]))
graph6.add(right(["o1", "o2"]))
graph6.add(fermion(["i1", "v1", "v2", "i2"]))
graph6.add(fermion(["o2", "v3", "v4", "o1"]))
graph6.add(photon(["v1", "v3"]))
graph6.add(photon(["v2", "v4"]))

test.write("./scratch/test.tex")
print(test)

test.compile("./scratch/test.tex", "test_mp.mp", True)
