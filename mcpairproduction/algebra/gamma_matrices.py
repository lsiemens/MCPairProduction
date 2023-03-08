"""Gamma Matrices

"""

import numpy

class cmatrix(numpy.ndarray):
    def __str__(self):
        if len(self.shape) > 2:
            print("high rank tensor")
            return super().__str__()
        if len(self.shape) < 2:
            obj = numpy.array([self]).view(type(self))
            return str(obj)

        # use the object dtype not str_ to get an array of arbitrary
        # length strings.
        elements = numpy.empty(shape=self.shape, dtype=object)

        n, m = self.shape
        width = 2

        for i in range(n):
            for j in range(m):
                elements[i, j] = f"{self[i, j]:.1g}"

#        mapped = numpy.zeros(shape=self.shape)
#        for key, value in mappings.items():
#            map = numpy.isclose(self, value)
#            mapped = numpy.logical_or(map, mapped)
#            elements[map] = " "*(width - len(key)) + key

#        unmapped = numpy.logical_not(mapped)
#        elements[unmapped] = numpy.rint(self[unmapped]).astype(str)

        

        line_width = self.shape[1]*(width + 1) - 1
        string = "┌" + " "*line_width + "┐\n"
        for row in elements:
            string += "│" + ",".join(row) + "│\n"
        string += "└" + " "*line_width + "┘"
        return string

class Repersentation:
    def __init__(self, generators):
        self.generators = [numpy.array(generator).view(cmatrix) for generator in generators]

    def __getitem__(self, key):
        return self.generators[key]

    def __str__(self):
        generators = [str(generator) for generator in self.generators]
        string = "Group/Algebra repersentation with the generators:\n" + "\n".join(generators)
        return string

gamma = Repersentation([ [[  1,  0,  0,  0],
                          [  0,  1,  0,  0],
                          [  0,  0, -1,  0],
                          [  0,  0,  0, -1]],

                         [[  0,  0,  0,  1],
                          [  0,  0,  1,  0],
                          [  0, -1,  0,  0],
                          [ -1,  0,  0,  0]],

                         [[  0,  0,  0,-1j],
                          [  0,  0, 1j,  0],
                          [  0, 1j,  0,  0],
                          [-1j,  0,  0,  0]],

                         [[  0,  0,  1,  0],
                          [  0,  0,  0, -1],
                          [ -1,  0,  0,  0],
                          [  0,  1,  0,  0]] ])

print(gamma)

print("\n\n\n")

print("square")
for i in range(4):
    print(f"gamma_{i}^2")
    print(numpy.matmul(gamma[i], gamma[i] + 0.1))

"""
print("product")
for i in range(4):
    for j in range(4):
        print(f"gamma_{i} * gamma_{j}")
        print(numpy.matmul(gamma[i], gamma[j]))

print("commutator")
for i in range(4):
    for j in range(4):
        print(f"[gamma_{i}, gamma_{j}]")
        print(numpy.matmul(gamma[i], gamma[j]) - numpy.matmul(gamma[j], gamma[i]))

print("anticommutator")
for i in range(4):
    for j in range(4):
        print(f"{{gamma_{i}, gamma_{j}}}")
        print(numpy.matmul(gamma[i], gamma[j]) + numpy.matmul(gamma[j], gamma[i]))

"""
