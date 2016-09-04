'''
Copyright (c) 2016, Tian Jiang
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

# Working version of adapt search with approximate results. Corresponding to peptoid6.pdf and adapted for a general rigid k-polygon by matching all (k choose 2) pairs of vertices' distances. This is not efficient for equilateral triangles. Using backbone 07AA1-6-C.pdb.

import numpy as np
from numpy import linalg as LA
import scipy, scipy.spatial
import itertools
import time
from scipy.sparse import dok_matrix
from octtree_minCube import *
from adaptfuns_suff import *
from searchcycles_suff import *
from choose import *

def peptoidDesign(err, eta, n, nsample, ncycles):
	start_time = time.time()
	##########################Params#######################################33
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
	#kgon = 3 # the target is a k-polygon, for k >= 3.
	#sidelen =[4.00,4.05,4.10] # the target polygons's side lengths ordered by the pair of vertices in the order Choose(kgon,2). 
	kgon = 3
	sidelen = [4.4, 4.1, 5.]
	# if all side chain bonds move by 5 degrees, the binding nodes move by 0.0785
	assert(kgon*(kgon-1)/2 == len(sidelen))
	#nsample = 8# sample nsample times in the free rotation of each node. 
	#print 'nsample = ', nsample
	#ncycles = 1 # number of desirable 3-tuples to find 
	#print 'ncycles = ', ncycles
	mpts = 1 # number of points to consider from each desirable cube
	#print 'mpts = ', mpts
	#err = 0.1   # target error tolerance
	#err = 0.5
	#print 'err = ', err
	#eta = 0.2 # approximation margin of width eta * err
	#print 'eta = ', eta
	minDist = 1.0
	#minDist = 2.2 # minimum distance between two nodes either from side chains or the backbone
	#print 'minDist between any two side chains/backbone =', minDist
	draw = False  # whether to draw the desirable polygons
	#print 'draw the desirable polygon = ', draw
	output = False
	#print 'produce outputs = ', output
	minIonDist = 1. # minimum distance between the binding ion and any of the backbone atoms
	# The angle, orientation and position of the side chains at the start from the backbone, data based on 07AA1-6-C.pdb
	#TODO: read the pdb file and generate the following array automatically
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
	CA1 = [3.355,2.419,7.397]
	C1 = [3.826,3.250,6.200]
	O1 = [4.683,2.827,5.422]
	CA2 = [4.005,5.409,5.173]
	C2 = [3.524,5.377,3.745]
	O2 = [2.332,5.227,3.475]
	CA3 = [5.878,5.621,3.029]
	C3 = [6.372,6.855,3.776]
	O3 = [5.664,7.863,3.846]
	CA4 = [8.399,5.579,4.333]
	C4 = [7.962,4.630,5.451]
	O4 = [7.146,4.979,6.305]
	CA5 = [7.878,2.419,6.357]
	C5 = [8.339,2.583,7.793] 
	O5 = [9.534,2.692,8.051] 
	CA6 = [5.974,2.539,8.503]  
	C6 = [5.419,1.199,7.994] 
	O6 = [6.116,0.222,7.995] 
	Ns = [N1,N2,N3,N4,N5,N6] # Take-off N's  
	TOCs = [sC1,sC2,sC3,sC4,sC5,sC6]  # Take-Off Carbons from N
	BB1 = [N1,CA1,C1,O1] # residue 1
	BB2 = [N2,CA2,C2,O2] # residue 2
	BB3 = [N3,CA3,C3,O3] # residue 3
	BB4 = [N4,CA4,C4,O4] # residue 4
	BB5 = [N5,CA5,C5,O5] # residue 5
	BB6 = [N6,CA6,C6,O6] # residue 6
	backbone = BB1+BB2+BB3+BB4+BB5+BB6 #backbone atoms without take-off carbons
	with open("6mer.pdb") as backbonefile:
		backbonelines = backbonefile.readlines()
	#n = 3 # number of carbon atoms on chain N-C-C-...-C-S
	totlen = nsample**n #total number of points on the manifold of sulphur atom
	Phi = np.arange(0.0, 2*np.pi, 2*np.pi/nsample)
	nPhis = []
	for i in xrange(n):
		nPhis.append(Phi)
	loopPhi = list(itertools.product(*nPhis))	
	# the 3-dimensional arrays that stores all the points
	# first dim is side chain index, second dim is vector size, third dim is # of total sampled points
	ptMat = np.ndarray(shape=(len(Ns),3,totlen),dtype=float,order='F')
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
	manifolds = [[] for x in range(len(Ns))]
	for chain in range(len(Ns)):
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
	position = [0 for x in range(len(Ns))]
	size = [0 for x in range(len(Ns))]
	for chain in range(len(Ns)):
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
	for chain in range(len(Ns)):
		myTree = OctTreeBox(np.array(position[chain]), initHalfSize, minSize)
		for i in range(totlen):
			pointOb = Point(chain, i, ptMat[chain,:,i])
			myTree.AddPoint(pointOb)
		treesList.append(myTree)
	End = time.time() - Start
	#print 'Octrees generated in %.2f seconds' %(End)
	#===========================================================================
	# Function takes pairs and prune two octrees to contain only points that appear in pairs
	# Input : pairs is a list of 2-tuples (p,q)
	# Output : octrees nt1,nt2 that are pruned t1,t2
	# Global effect : None
	def PruneTrees(pairs, t1, t2):
		sc1 = set()
		sc2 = set()
		for (cp,cq) in pairs:
			sc1.add(cp)
			sc2.add(cq)
		nt1 = Prune(sc1, t1)
		nt2 = Prune(sc2, t2)
		return nt1,nt2

	# Function builds a new tree whose leaves are "leaves" and is a subtree of "tree"
	# Input : leaves is a set of leaf cubes that are contains in the octree tree
	# Output : a pruned ntree
	# global effect : none
	def Prune(leaves, tree):
		if tree.isLeaf:
			if tree in leaves:
				return tree
			else:
				return None
		else:
			newChildren = [Prune(leaves, child) for child in tree.GetChildren()]
			newChildren = [x for x in newChildren if x != None]
			if newChildren == []:
				return None
			else:
				newtree = OctTreeBox(tree.center, tree.halfSize, minSize)
				newtree.isLeaf = False
				newtree.children = newChildren
				newtree.points = None
				return newtree

	#===========================================================================
	scs = Choose(len(Ns),kgon) # a list of takeoff sites lists
	scs_order = []
	for i in range(len(scs)):
		scs_order += Arrange(set(scs[i]))
	lenPairs = Choose(kgon,2) # a list of length pairs in the same order of sidelen
	totalsearchstart = time.time()
	for inds in scs_order:
		#print "=========== Taking off from ", inds, " ==============="
		trees = []
		for i in range(len(inds)):
			trees.append(treesList[inds[i]])
		nPairs = []
		Search = True
		startm = time.time()
		for i in range(len(lenPairs)):
			[i1,i2] = lenPairs[i]
			'''
			# crude check for initial cubes against necessary condition
			initCs = LA.norm(np.array(position[inds[L[0]]]) - np.array(position[inds[L[1]]]))
			initSize = initHalfSize * 2
			if (initCs - (3.**0.5)*initSize > sidelen[i]) or (initCs + (3.**0.5)*initSize < sidelen[i]):
				print 'Initial check screened out this possibility!'
				Search = False
				break
			'''
			adaptStartm = time.time()
			pairs = AdaptSearchApprox(trees[i1],trees[i2], sidelen[i], err)
			#print [inds[i1],inds[i2]], ' : found %d pairs, used %.2f seconds' %(len(pairs), time.time() - adaptStartm)
			sys.stdout.flush()
			if len(pairs) == 0:
				Search = False
				break
			nPairs.append(pairs)
			trees[i1],trees[i2] = PruneTrees(pairs, trees[i1], trees[i2])
			#print pairs[0][0].isLeaf, pairs[0][0].halfSize*2*(3**.5)
		if Search:
			# nPairs = Trim(nPairs)
			cycleStart = time.time()
			cycles = SearchPolys(nPairs, kgon, manifolds, nsample, mpts, minDist, minIonDist, Ns, backbone, backbonelines, ncycles, draw, output)
			#print "Total time spent finding all %d desirable k-gons is %.2f seconds." %(len(cycles), time.time() - cycleStart)
			#print "There are no desirable k-gons since at least one desirable pair of vertices is missing."
		#endm = time.time() - startm
		#print 'Total time spent for this target polygon is %.2f seconds' %(endm)
	print 'err=%.2f, eta=%.1f, n=%d, nsample=%d, ncycles=%d:' %(err, eta, n, nsample, ncycles)
	print 'Total time: %.2f' %(time.time() - start_time)

if __name__ == '__main__':
	for err in [0.05, 0.1, 0.5]:
		for eta in [0.1, 0.2, 0.3]:
			for n in [3,2]:
                                peptoidDesign(err, eta, n, 8, 1)
