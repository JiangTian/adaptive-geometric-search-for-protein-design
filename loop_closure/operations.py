'''
Copyright (c) 2016, Tian Jiang
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import numpy as np

jokerElem = -1

# long term TODO: all negative values except -1 are wasted.
# encode angle and x,y,z into a unique integer
def code(angle, x, y, z, ba, bb):
    return angle + ba * (x + bb * (y + bb * z)) 

# computes the angle from code
def angle(code, ba):
    return code % ba

# computes the x,y,z from code
def xyz(code, ba, bb):
    z = code // (ba * bb**2)
    r = code - z * (ba * bb**2) 
    y = r // (ba * bb)
    x = (r - y * (ba * bb)) // ba
    return x, y, z

# develop the semi-loop in "nsteps" many residues and considue only the branches within the designated space "pocket"
def forward(M, code0, nsteps, pocket, ba, bb):
    print "Forward pass. Ensure you're space is big enough, or it will produce incorrect results"
    assert(M.shape[0] == ba)
    nmaxpointsperangle = M.shape[1]
    C = [np.ndarray([1], dtype='int64')]
    C[0][0] = code0
    for i in range(1, nsteps):
        print "Step %d, %d starting points."%(i, C[-1].shape[0])
        npoints = C[i-1].shape[0]
        C2 = np.ndarray([npoints, nmaxpointsperangle], dtype='int64')
        for j in range(npoints):
            c = C[i-1][j]
            a = angle(c, ba)
            np.add(M[a], c - a, C2[j])
            isJoker = (M[a] == jokerElem) #TODO can be better...
            np.multiply(C2[j], 1-isJoker, C2[j])
            C2[j] = C2[j] - isJoker
        #TODO: check if we get out of space
        C.append(np.unique(C2))
        if C[-1][0] == -1:
            # remove joker
            C[-1] = C[-1][1:]
        C[-1] = prune(C[-1], pocket, ba, bb)
    print "Found %d 5d points"%(C[-1].shape[0])
    return C

# search for paths backward to see results
# TODO: precompute a reverse lookup table of M to accelerate production of many results
def backward(M, C, m, nsteps, ba):
    print "Backward pass. Produce a result."
    chain = [-1 for i in range(nsteps)]
    chain[nsteps-1] = m
    for i in range(nsteps-2, -1, -1):
        p = chain[i+1]
        print "Step %d, %d starting points."%(i, C[i].shape[0])
        npoints = C[i].shape[0]
        for j in range(npoints):
            c = C[i][j]
            a = angle(c, ba)
            nextBoxes = M[a] + c - a
            if p in nextBoxes:
                chain[i] = c
                break
        assert(chain[i] != -1)
    return np.array(chain)

# delete all those branches that fall out of "pocket"
def prune(l, pocket, ba, bb):
    x, y, z = xyz(l, ba, bb)

    xbool = (pocket[0][0] <= x) * (x <= pocket[0][1])
    ybool = (pocket[1][0] <= y) * (y <= pocket[1][1])
    zbool = (pocket[2][0] <= z) * (z <= pocket[2][1])

    allbool = xbool * ybool * zbool
    return l[allbool]

# match cubes and angles in the middle where two semi-loops meet
def matching(c1, c2, match, nSample, ba):
    print "Matching... %d vs %d 5d points"%(c1.shape[0], c2.shape[0])
    # 1) match positions
    print "Matching positions"
    angles1 = angle(c1, ba)
    boxes1 = c1 - angles1
    angles2 = angle(c2, ba)
    boxes2 = c2 - angles2
    intersection = set(np.unique(boxes1)).intersection(set(np.unique(boxes2)))
    print (len(intersection))
    # 2) match angles
    print "Matching angles"
    m1, m2 = [], []
    anglesPerBox1 = {x:[] for x in intersection}
    anglesPerBox2 = {x:[] for x in intersection}
    for i in range(angles1.shape[0]):
        if boxes1[i] in intersection:
            anglesPerBox1[boxes1[i]].append(angles1[i])
    for i in range(angles2.shape[0]):
        if boxes2[i] in intersection:
            anglesPerBox2[boxes2[i]].append(angles2[i])
    for box in intersection:
        ang1 = anglesPerBox1[box]
        ang2 = anglesPerBox2[box]
        for a1 in ang1:
            ang2_match = np.array(match[a1 * nSample]) // nSample
            A2 = set(np.unique(ang2_match)).intersection(set(ang2))
            m1 += [a1 + box for i in range(len(A2))]
            m2 += [a2 + box for a2 in A2]
            return np.array(m1), np.array(m2)

# compute and return M, the lookup table for matching residue conformations 
def getM(conformations, match, nSample, nStart, startbox, boxsize, ba, bb): 
    M2 = [[] for i in range(nStart)]
    for ind1 in range(nSample * nStart):
        inAngle = ind1 // nSample
        for ind2 in match[ind1]:
            outAngle = ind2 // nSample
            outPos = list((conformations[0][-1, ind2, :] + boxsize/2.) // boxsize)
            outPos = map(int, outPos)
            M2[inAngle].append([outAngle] + outPos)
    M3 = [np.unique(np.array([code(m[0], m[1], m[2], m[3], ba, bb) for m in x], dtype='int64')) for x in M2]
    nmaxpointsperangle = max([len(x) for x in M3])
    
    M = np.ndarray([ba, nmaxpointsperangle], dtype = 'int64')
    M.fill(-1)
    for i in range(len(M3)):
        m = M3[i]
        for j in range(len(m)):
            M[i][j] = m[j]
    return M
