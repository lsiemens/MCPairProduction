"""Gamma Matrices

"""

import numpy

class gamma:
    def __init__(self):
        self._gamma = numpy.array([ [[  1,  0,  0,  0],
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
                                     [  0,  1,  0,  0]] ]

print(gamma)

for mat in gamma:
    print(numpy.array_str(numpy.matmul(mat, mat), suppress_small=True))
