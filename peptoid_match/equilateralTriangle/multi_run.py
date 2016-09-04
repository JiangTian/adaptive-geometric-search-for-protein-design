'''
Copyright (c) 2016, Tian Jiang
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

# Adaptive search with target polygon being an equilateral triangle.

import numpy as np
from numpy import linalg as LA
import scipy, scipy.spatial
import itertools
import time
from scipy.sparse import dok_matrix
from octtree_minCube import *
from adaptfuns_suff import *
from searchcycles_suff import *


def peptoidDesign(err, eta, n, nsample, ncycles):
	start_time = time.time()
	# The following values are the best bond lengths and angles for each Carbon-Carbon link or Carbon-Sulfur link from molecular links table
	thetaNCC = np.pi*(180.0-110.0)/180.0 # pi minus the best NCC angle
	bondNC = 1.455
	thetaCCC = np.pi*(180.0-120.0)/180.0
	bondCC = 1.5
	thetaCCS = np.pi*(180.0-114.5)/180.0
	bondCS = 1.817
	bonds = [bondCC,bondCS]
	lengths = [bondCC] * (n-1) + [bondCS]
	thetas = [thetaNCC] + [thetaCCC] * (n-2) + [thetaCCS]
	sidelen = 4.05 # the desired equilateral triangle's side length
	#nsample = 10 # sample nsample times in the free rotation of each node. 
	#print 'nsample = ', nsample
	#ncycles = 1 # number of desirable 3-tuples to find 
	draw = False
	output = False
	#err = sidelen*0.01    # target error tolerance 1%
	#eta = 0.2 # approximation margin of width eta * err
	minDist = 1. # minimum distance between two nodes either from side chains or the backbone
	minIonDist = 2.5 # minimum distance between the binding ion and any of the backbone atoms
	# The angle, orientation and position of the side chains at the start from the backbone, data based on 07AA1-6-C.pdb
	# TODO: pipeline data reading 
	N1 = [4.134,1.199,7.560]
	sC1 = [3.558,-0.033,7.026]
	N2 = [3.300,4.494,6.068]
	sC2 = [2.283,5.065,6.951]
	N3 = [4.441,5.645,2.777]
	sC3 = [3.955,5.949,1.434]
	N4 = [7.598,6.793,4.339]
	sC4 = [8.053,7.901,5.174]
	N5 = [8.485,3.385,5.442]
	sC5 = [9.437,2.924,4.429]
	N6 = [7.399,2.480,8.768]
	sC6 = [7.851,2.166,10.120]
	Ns = [N1,N2,N3,N4,N5,N6]  
	TOCs = [sC1,sC2,sC3,sC4,sC5,sC6]  # Take-Off Carbons
	with open("6mer.pdb") as backbonefile:
		backbonelines = backbonefile.readlines()
	#n = 3 # number of carbon atoms on chain N-C-C-...-C-S
	totlen = nsample**n #total number of points on the manifold of sulphur atom
	Phi = np.arange(0.0, 2*np.pi, 2*np.pi/nsample)
	nPhis = []
	for i in xrange(n):
		nPhis.append(Phi)
	loopPhi = list(itertools.product(*nPhis))	
	# the 3-dimensional arrays that stores all the points and their phi positions by layers
	# first dim is side chain index, second dim is vector size, third dim is # of total sampled points
	ptMat = np.ndarray(shape=(6,3,totlen),dtype=float,order='F')
	posMat = np.ndarray(shape=(6,n,totlen),dtype=float,order='F')

	########################### Custom Functions ############################################

	# root function computes the spherical angle coordinates of the N-C link, first link on side chain
	# Input : cartesian coordinates of N and C. C is the take off(first) carbon on the side chain.
	# Output : theta and phi, spherical coordinate of the first C node on the side chain
	# Global side effects : None
	def root(N,C):
		V = np.array(C)-np.array(N)
		theta = np.arccos(V[2]/LA.norm(V))
		if V[1]<0 and V[0]>0:
			phi = 2*np.pi - np.arctan(-V[1]/V[0])
		elif V[1]<0 and V[0]<0:
			phi = np.pi + np.arctan(V[1]/V[0])
		elif V[1]>0 and V[0]<0:
			phi = np.pi - np.arctan(-V[1]/V[0])
		else:
			phi = np.arctan(V[1]/V[0])
		return [theta,phi]
	# calculate the initial tilt of each side chain
	Ch1 = root(Ns[0],TOCs[0])
	Ch2 = root(Ns[1],TOCs[1])
	Ch3 = root(Ns[2],TOCs[2])
	Ch4 = root(Ns[3],TOCs[3])
	Ch5 = root(Ns[4],TOCs[4])
	Ch6 = root(Ns[5],TOCs[5])
	root_theta = [Ch1[0], Ch2[0], Ch3[0], Ch4[0], Ch5[0], Ch6[0]] 
	root_phi = [Ch1[1], Ch2[1], Ch3[1], Ch4[1], Ch5[1], Ch6[1]] 

	#=======================================================================
	#print '...generating manifolds using matrix multiplication'
	def rotateM(phi, theta):
		u = (np.cos(phi+np.pi/2.0),np.sin(phi+np.pi/2.0),0.0)
		rotM = np.array([[np.cos(theta)+u[0]**2*(1-np.cos(theta)), u[0]*u[1]*(1-np.cos(theta)), u[1]*np.sin(theta)],
		[u[0]*u[1]*(1-np.cos(theta)), np.cos(theta)+u[1]**2*(1-np.cos(theta)), -u[0]*np.sin(theta)],
		[-u[1]*np.sin(theta), u[0]*np.sin(theta), np.cos(theta)]])
		return rotM

	def GetRotations(l, theta, nsample):
		v = np.zeros([nsample, 3])
		R = np.zeros([nsample, 3, 3])
		d = l*np.sin(theta)
		z = l*np.cos(theta)
		#TODO: store sin, cos values to optimize
		for i in range(nsample):
			phi = Phi[i]
			v[i] = [d * np.cos(phi), d * np.sin(phi), z]
			R[i] = rotateM(phi, theta)
		return R, np.matrix(v)

	# R0 = rotation matrix of the NC
	# v0 = position of the take-off C
	def BuildManifold(lengths, thetas, nsample, R0, v0):
		csLength = len(thetas)
		manifolds = [None for x in range(csLength + 1)]
		manifolds[0] = np.matrix(np.zeros([3, 1]))
		for i in range(csLength):
			R, v = GetRotations(lengths[csLength-i-1], thetas[csLength-i-1], nsample)
			for j in range(i,-1,-1):
				M = manifolds[j]
				m = M.shape[1]
				M0 = np.matrix(np.zeros([3, m * nsample]))
				#TODO: can be one big matrix multi with scrambling
				for k in range(nsample):
					M0[:, (k*m) : ((k+1)*m)] = np.matrix(R[k]) * M + np.transpose(v[k])
				manifolds[j+1] = M0
		for k in range(csLength + 1):
			manifolds[k] = np.matrix(R0) * manifolds[k] + np.transpose(v0)
		return manifolds

	startsamp = time.time()
	manifolds = [[] for x in range(6)]
	for chain in range(6):
		R0 = rotateM(root_phi[chain], root_theta[chain])
		v0 = np.matrix(TOCs[chain])
		manifolds[chain] = BuildManifold(lengths, thetas, nsample, R0, v0)
		ptMat[chain] = np.matrix(manifolds[chain][-1])
	endsamp = time.time() - startsamp
	#print 'Time spent generating manifolds is %.2f seconds' %(endsamp)

	#============= generate Octrees ==========================================
	#print '...generating octrees for all manifolds'
	Start = time.time()
	treesList = []
	position = [0,0,0,0,0,0]
	size = [0,0,0,0,0,0]
	for chain in range(6):
		xmin = min(ptMat[chain,0,:])
		xmax = max(ptMat[chain,0,:])
		ymin = min(ptMat[chain,1,:])
		ymax = max(ptMat[chain,1,:])
		zmin = min(ptMat[chain,2,:])
		zmax = max(ptMat[chain,2,:])
		position[chain] = [(xmin+xmax)/2.0, (ymin+ymax)/2.0, (zmin+zmax)/2.0]
		size[chain] = max(xmax-xmin, ymax-ymin, zmax-zmin)
	initHalfSize = max(size)/2.0
	minSize = eta*err/(4.0*(3.**0.5))
	for chain in range(6):
		myTree = OctTreeBox(np.array(position[chain]), initHalfSize, minSize)
		for i in range(totlen):
			pointOb = Point(chain, i, ptMat[chain,:,i])
			myTree.AddPoint(pointOb)
		treesList.append(myTree)
	End = time.time() - Start
	#print 'Octrees generated in %.2f seconds' %(End)

	#===========================================================================
	subsets = [[0,1,2],[0,1,3],[0,1,4],[0,1,5],[0,2,3],[0,2,4],[0,2,5],[0,3,4],[0,3,5],[0,4,5],[1,2,3],[1,2,4],[1,2,5],[1,3,4],[1,3,5],[1,4,5],[2,3,4],[2,3,5],[2,4,5],[3,4,5]]

	avepts = [] # list of np arrays
	for inds in subsets:
		#print "=========== Taking off from N %d, %d, %d ===========" %(inds[0],inds[1],inds[2])
		nPairs = []
		Prod = 1
		startm = time.time()
		for i in range(3):
			P = dok_matrix((totlen, totlen), dtype=np.float32)
			adaptStartm = time.time()
			pairs, numpts = AdaptSearchApprox(treesList[inds[i]],treesList[inds[(i+1)%3]], sidelen, err)
			avepts.append(np.array(numpts))
			#print 'AdaptSearch found %d pairs, used %.2f seconds' %(len(pairs), time.time() - adaptStartm)
			nPairs.append(pairs)
		cycleStart = time.time()
		cycles = SearchCyclesApprox(nPairs, manifolds, nsample, minDist, minIonDist, Ns, backbonelines, ncycles, draw, output)
		#print "Total time spent finding all %d cycles is %.2f seconds." %(len(cycles), time.time() - cycleStart)
	end_time = time.time() - start_time
	print 'err=%.2f, eta=%.1f, n=%d, nsample=%d, ncycles=%d:' %(err, eta, n, nsample, ncycles)
	print 'Total time: %.2f' %(end_time)

if __name__=='__main__':
	for err in [0.05, 0.1, 0.5]:
		for eta in [0.1, 0.2, 0.3]:
			for n in [3, 2]:	
				peptoidDesign(err, eta, n, 8, 1)
