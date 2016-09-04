'''
Copyright (c) 2016, Tian Jiang
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import numpy as np 
import numpy.linalg as LA
import math, random
import copy

# produce nearly uniformly distanced points on a sphere
def fibonacci_sphere(samples, randomize=False):
    rnd = 1.
    if randomize:
        rnd = random.random() * samples

    points = []
    offset = 2./samples
    increment = math.pi * (3. - math.sqrt(5.))

    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2)
        r = math.sqrt(1 - pow(y,2))
        phi = ((i + rnd) % samples) * increment
        x = math.cos(phi) * r
        z = math.sin(phi) * r
        points.append([x,y,z])

    return points

def getAngle(u,v):
    return np.arccos(np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)))

# hardcoded for now
def getBondStat(side): #TODO: data processing
    N = np.array([95.709,  34.056,  35.641])
    CA = np.array([94.301,  34.348,  35.958])
    C = np.array([93.914 , 34.453,  35.727])
    Nl = np.array([ 93.619,  34.376,  35.587])
    CAl = np.array([93.414,  34.389,  35.264])
    angles = [] # Ca, C, N, Ca_next, ..
    if side == 1:
        angles.append(getAngle(CA-C, Nl-C))
        angles.append(getAngle(C-Nl, CAl-Nl))
        angles.append(getAngle(N-CA, C-CA))
    else:
        angles.append(getAngle(CAl-Nl, C-Nl))
        angles.append(getAngle(Nl-C, CA-C))
        angles.append(getAngle(N-CA, C-CA))
    lengths = []
    if side == 1:
        lengths.append(np.linalg.norm(C-CA))
        lengths.append(np.linalg.norm(Nl-C))
        lengths.append(np.linalg.norm(CAl-Nl))
    else:
        lengths.append(np.linalg.norm(CAl-Nl))
        lengths.append(np.linalg.norm(Nl-C))
        lengths.append(np.linalg.norm(C-CA))
    return angles, lengths

# rotation matrix from (1,0,0) to pt
def rotate(pt):
    rot = np.array([[pt[0], 0., pt[1]**2 + pt[2]**2], 
                    [pt[1], -pt[2], -pt[0] * pt[1]],
                    [pt[2], pt[1], -pt[0] * pt[2]]])
    return rot

def get3DTensor(res, nStart, nSample, side):
    uniformPts = fibonacci_sphere(nStart)
    angles, lengths = getBondStat(side)
    assert(len(lengths) == 3) # Doug said all residues have length 3
    C = lengths[0] * np.array(uniformPts)
    layer1 = np.repeat(C, nSample, axis = 0)
    
    pts = np.array([[lengths[1]*math.cos(math.pi-angles[0]),
                     lengths[1]*math.sin(math.pi-angles[0])*math.sin(2.*math.pi*ibeta/nSample),
                     lengths[1]*math.sin(math.pi-angles[0])*math.cos(2.*math.pi*ibeta/nSample)]
                    for ibeta in range(nSample)]) # sample points around (1,0,0) with angle angles[0]
    pts += np.array([lengths[0], 0, 0])
    rotations = map(rotate, uniformPts)
    
    # generate Ns
    layer2 = []
    for rot in rotations:
        layer2 += [np.dot(rot, pt) for pt in pts]
    layer2 = np.array(layer2)

    # Assume first and last bond angles are the same, then Ca_next = N + scaled(C)
    layer3 = layer1*(lengths[2]/lengths[0]) + layer2 
    
    tensor = np.ndarray([len(lengths), layer1.shape[0], 3])
    tensor[0] = layer1
    tensor[1] = layer2
    tensor[2] = layer3
    return tensor
    
# stores all residue conformations for each residue ID
def buildConformations(resIDs, nStart, nSample, side):    
    return {res : get3DTensor(res, nStart, nSample, side) for res in resIDs}

# during semi-loop development, returns a vector of all matching residue indices for residue index i
def checkAngle(i, totalInd, targetAngle, angleErr, conformations, resID):
    u = conformations[resID][-2, i, :] - conformations[resID][-1, i, :]
    Vs = conformations[resID][0, :, :]
    L = LA.norm(u) * LA.norm(Vs, axis = 1)
    cosVal = np.clip(np.dot(Vs, u) / L, -1, 1)
    angles = np.arccos(cosVal)
    boolvec = abs(angles - targetAngle) < angleErr
    return [i for i, x in enumerate(boolvec) if x == 1] #TODO: make it faster

# returns a vector of all matching residue indices for residue index i in the middle of the loop
def checkAngle_mid(i, totalInd, targetAngle, angleErr, conformations1, conformations2, resID):
    u = conformations1[resID][-2, i, :] - conformations1[resID][-1, i, :]
    Vs = conformations1[resID][-2, :, :] - conformations1[resID][-1, :, :]
    L = LA.norm(u) * LA.norm(Vs, axis = 1)
    cosVal = np.clip(np.dot(Vs, u) / L, -1, 1)
    angles = np.arccos(cosVal)
    boolvec = abs(angles - targetAngle) < angleErr
    return [i for i, x in enumerate(boolvec) if x == 1] 

def plotCube(ax, pocket):
    ax.plot([pocket[0][0], pocket[0][1]], [pocket[1][0], pocket[1][0]], [pocket[2][0], pocket[2][0]], c='black')
    ax.plot([pocket[0][0], pocket[0][1]], [pocket[1][0], pocket[1][0]], [pocket[2][1], pocket[2][1]], c='black')
    ax.plot([pocket[0][0], pocket[0][1]], [pocket[1][1], pocket[1][1]], [pocket[2][0], pocket[2][0]], c='black')
    ax.plot([pocket[0][0], pocket[0][1]], [pocket[1][1], pocket[1][1]], [pocket[2][1], pocket[2][1]], c='black')
    ax.plot([pocket[0][0], pocket[0][0]], [pocket[1][0], pocket[1][1]], [pocket[2][0], pocket[2][0]], c='black')
    ax.plot([pocket[0][0], pocket[0][0]], [pocket[1][0], pocket[1][1]], [pocket[2][1], pocket[2][1]], c='black')
    ax.plot([pocket[0][1], pocket[0][1]], [pocket[1][0], pocket[1][1]], [pocket[2][0], pocket[2][0]], c='black')
    ax.plot([pocket[0][1], pocket[0][1]], [pocket[1][0], pocket[1][1]], [pocket[2][1], pocket[2][1]], c='black')
    ax.plot([pocket[0][0], pocket[0][0]], [pocket[1][0], pocket[1][0]], [pocket[2][0], pocket[2][1]], c='black')
    ax.plot([pocket[0][0], pocket[0][0]], [pocket[1][1], pocket[1][1]], [pocket[2][0], pocket[2][1]], c='black')
    ax.plot([pocket[0][1], pocket[0][1]], [pocket[1][0], pocket[1][0]], [pocket[2][0], pocket[2][1]], c='black')
    ax.plot([pocket[0][1], pocket[0][1]], [pocket[1][1], pocket[1][1]], [pocket[2][0], pocket[2][1]], c='black')


    
