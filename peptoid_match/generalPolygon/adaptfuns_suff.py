'''
Copyright (c) 2016, Tian Jiang
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

from octtree_minCube import *
import numpy as np
import math
from itertools import product


def norm2(x):
    return np.dot(x, x)

# returns a list of 2-tuples
def AdaptSearchApprox(t1, t2, l, err):
    result = []
    assert(t1.GetDepth() == t2.GetDepth())
    depth = t1.GetDepth()
    combos = [[] for x in range(depth+1)]
    combos[0] = [(t1, t2)]
    for i in range(depth):
        for (b0, b1) in combos[i]:
            combos[i+1] += Compare_nece(b0, b1, l, err)
    for (b0, b1) in combos[depth]:
	result += Compare_suff(b0, b1, l, err)
    return result

def Compare_nece(cube0, cube1, l, err):
    # cube0.halfSize = 2 * boxi.halfSize
    threshold = err + math.sqrt(3)*cube0.halfSize
    tau0 = (max(0., l - threshold)) ** 2
    tau1 = (l + threshold) ** 2
    def isGood(pi, pj):
        n = norm2(pi.center - pj.center)
        return (tau0 <= n) and (n <= tau1)
    return [(boxi, boxj) for (boxi, boxj) in product(cube0.GetChildren(), cube1.GetChildren()) if isGood(boxi, boxj)]

def Compare_suff(cube0, cube1, l, err):
    threshold = err - math.sqrt(3)*cube0.halfSize*2
    tau0 = (max(0., l - threshold)) ** 2
    tau1 = (l + threshold) ** 2
    n = norm2(cube0.center - cube1.center)
    if ((tau0 <= n) and (n <= tau1)):
    	return [(cube0, cube1)]
    else:
	return []
