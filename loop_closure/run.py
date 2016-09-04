'''
Copyright (c) 2016, Tian Jiang
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import numpy as np
import support
import operations
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

starttime = time.time()
# set params:
loop_length = 12
boxsize = 0.1
pivotShift = np.array([5.0, 0.0, 0.0])
start = [0., 0., 0.]
range3D = [121, 121, 121]
startbox = np.array(range3D)/2
currentLeftLayer = np.array([startbox])
resIDs = [0] * (loop_length/2) # TODO : for now all resIDs are the same
nStart = 100
nSample = 10 # TODO : pick nStart, nSample and boxsize to bound errors
totalInd = nStart * nSample
ba = nStart # encoding of angle
bb = max(range3D) # encoding of xyz
angles1, lengths1 = support.getBondStat(1)
angles2, lengths2 = support.getBondStat(2)
targetAngle1 = angles1[-1]
targetAngle2 = angles2[-1]
angleErr = 0.2
pocket = [[range3D[0]/2 - 50, range3D[0]/2 + 5],
          [range3D[1]/2 - 50, range3D[1]/2 + 5], 
          [range3D[2]/2 - 50, range3D[2]/2 + 5]]

print "1. Building conformations.."
conformations1 = support.buildConformations(resIDs, nStart, nSample, 1)
conformations2 = support.buildConformations(resIDs, nStart, nSample, 2)

print ".. and presearch good conformations combos."
match1 = {i:[] for i in range(totalInd)}
for i in range(totalInd):
        match1[i] = support.checkAngle(i, totalInd, targetAngle1, angleErr, conformations1, 0) #resID = 0
match2 = {i:[] for i in range(totalInd)}
for i in range(totalInd):
        match2[i] = support.checkAngle(i, totalInd, targetAngle2, angleErr, conformations2, 0) #resID = 0
match_mid = {i:[] for i in range(totalInd)}
for i in range(totalInd):
        match_mid[i] = support.checkAngle_mid(i, totalInd, targetAngle1, angleErr, conformations1, conformations2, 0) #resID = 0
print "matches done."

print "2. getting M"
M1 = operations.getM(conformations1, match1, nSample, nStart, startbox, boxsize, ba, bb)
M2 = operations.getM(conformations2, match2, nSample, nStart, startbox, boxsize, ba, bb)

print "checking the space is big enough"
#needed_space = (loop_length / 2 + 1) * sum(lengths) / boxsize
needed_space = (loop_length / 2 + 1) * sum(lengths1) / boxsize + 10#TODO: this is for debug
assert(needed_space <= range3D[0] and needed_space <= range3D[1] and needed_space <= range3D[2])
print "There is a TODO here."
#TODO: careful: here we assume the starting point is 0,0,0

start0 = np.array([range3D[0]/2, range3D[1]/2, range3D[2]/2])
start1 = start0 + pivotShift

growtime = time.time()
code0 = operations.code(0, start0[0], start0[1], start0[2], ba, bb) #TODO: origin angle
C1 = operations.forward(M1, code0, loop_length / 2 + 1, pocket, ba, bb) # TODO the -1 are fishy
print "%.2f seconds in growing one side" %(time.time() - growtime)

code0 = operations.code(1, start1[0], start1[1], start1[2], ba, bb) #TODO: origin angle
if loop_length % 2 == 0: 
	semiloop = loop_length / 2 + 1
else:
	semiloop = loop_length / 2 + 2
C2 = operations.forward(M2, code0, semiloop, pocket, ba, bb) # TODO the -1 are fishy
'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plotCube(ax, pocket)
ax.scatter([start0[0]], [start0[1]], [start0[2]], c='blue')
ax.scatter([start1[0]], [start1[1]], [start1[2]], c='green')
x, y, z = operations.xyz(C1[-1])
x, y, z = x[::100], y[::100], z[::100]
ax.scatter(x, y, z, c='purple')
x, y, z = operations.xyz(C2[-1])
x, y, z = x[::100], y[::100], z[::100]
ax.scatter(x, y, z, c='red')
plt.show()
'''
# N - CA - C, change the order of the other half of loop 
# check matching after developing N 
matchtime = time.time()
m1, m2 = operations.matching(C1[-1], C2[-1], match_mid, nSample, ba)
print "%.2f seconds in matching" %(time.time() - matchtime)

backtime = time.time()
for i in range(len(m1)):
    assert(m1[i] in C1[-1])
    chain1 = operations.backward(M1, C1, m1[i], loop_length/2 + 1, ba)

    assert(m2[i] in C2[-1])
    chain2 = operations.backward(M2, C2, m2[i], semiloop, ba)
    print "Total time used till here: %.2f seconds." %(time.time() - starttime)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plotCube(ax, pocket)
    ax.scatter([start0[0]], [start0[1]], [start0[2]], c='blue')
    ax.scatter([start1[0]], [start1[1]], [start1[2]], c='green')
    x1, y1, z1 = operations.xyz(chain1, ba, bb)

    '''
    angles1 = operations.angle(chain1)
    x_chain, y_chain, z_chain = [],[],[]
    for i in range(len(chain1)-1):
        box_shift = (np.array(conformations1[0][:, angles1[i+1]*nSample, :]) + boxsize/2.) // boxsize
        x_chain += list(x1[i] + np.array(box_shift[:, 0]))
        y_chain += list(y1[i] + np.array(box_shift[:, 1]))
        z_chain += list(z1[i] + np.array(box_shift[:, 2]))
    print x_chain
    ax.plot(x_chain, y_chain, z_chain, c='blue')
    '''

    x2, y2, z2 = operations.xyz(chain2, ba, bb)
    ax.plot(x1, y1, z1, c='blue')
    ax.plot(x2, y2, z2, c='green')
    plt.show()

print "Total time used: %.2f seconds." %(time.time() - starttime)

