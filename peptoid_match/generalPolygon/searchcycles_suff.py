'''
Copyright (c) 2016, Tian Jiang
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import numpy as np
import scipy, scipy.spatial
from numpy import linalg as LA
from octtree_minCube import *
import itertools
from output import *
from choose import *

def Graph(pairs):
	graph = {}
	for (i,j) in pairs:
		if i not in graph:
			graph[i] = []
		graph[i] += [j]
	return graph

# reverse ordered
def GraphR(pairs):
	graph = {}
	for (i,j) in pairs:
		if j not in graph:
			graph[j] = []
		graph[j] += [i]
	return graph

# function trims away all cubes that are loose ends iteratively
# Input : nPairs a list of n list of pair tuples
# Output : the number of elements removed or -1 if all were removed
# global effect : Trimmed nPairs[0]
def Trim(nPairs, kgon):
	if kgon > 1:
		n = kgon
		result = Trim(nPairs[-(len(nPairs)-(n-1)):], kgon-1)
		if result == 0:
			return 0
		sets = [set() for x in range(n-1)]
		for node in range(n-1):
			for i in range(len(nPairs[node])):
				sets[node].add(nPairs[node][i][0])
		valid = set.intersection(*sets)
		if len(valid) == 0:
			return 0
		for node in range(n-1):
			nPairs[node] = [p for p in nPairs[node] if p[0] in valid]
	else:
		return 1


# The main function in SearchPolys, it searches for all the desirable k-tuple cubes. 
# Input : nPairs gives all the pairwise desirable cubes in the order of Choose(len(Ns),kgon)
# Output : a list of k-list of cubes
# Global effect : None
def SearchGraph(pairs, k):	
	n = len(pairs)
	if k == 2:
		assert(n==1)
		return list(pairs[-1])
	else:
		result = []
		gs = [[] for x in range(k-1)]
		for i in range(1,k-1):
			gs[i] = Graph(pairs[i])
		gs[0] = GraphR(pairs[0])
		leng = n - (k - 1) 
		Cpost = SearchGraph(pairs[-leng:] , k-1) # gives all cube lists that 
		for cs in Cpost:
			if (cs[0] in gs[0]):
				for c0 in gs[0][cs[0]]:
					for i in range(k-2):
						if c0 in gs[i+1]:
							if cs[i+1] in gs[i+1][c0]:
								if i == k-3:
									result.append([c0]+list(cs))
							else: 
								break
						else:
							break
		return result

def PtsInCube(c,m):
	return c.GetPoints()[:m]

# Function SearchPolys loops through all possibilities and produces desirable k-gons
# Input : Npairs is a list of kgon lists of 2 tuples.
# e.g. Npairs = [[(1,2),(2,3)],[(0,1),(1,3)],[(2,3)]] 
# Output : a list of polygons
# Global effect : None
def SearchPolys(nPairs, kgon, manifolds, nsample, mpts, minDist, minIonDist, Ns, backbone, backbonelines, ncycles=None, draw=None, output=None):
	result = []
	CubeList = SearchGraph(nPairs, kgon)
	for kcube in CubeList:
		ptslist = [PtsInCube(cube,mpts) for cube in kcube]
		kpts = list(itertools.product(*ptslist))
		for kpt in kpts:
			assert(len(kpt) == kgon)
			nodes = [[] for x in range(kgon)]
			for i in range(kgon):
				nodes[i] = SideChainNodes(kpt[i], manifolds, nsample)
			if NoClash(nodes, Ns, backbone, minDist, minIonDist) and CheckAngle(nodes, np.pi/3):
				result.append(kpt)
				if draw == True:
					Draw(kpt,nodes,Ns)
				if output == True:
					Output(kpt,nodes,backbonelines,len(result))
				if ncycles!= None and len(result) == ncycles:
					return result
	return result


# Function TraceSideChain produces the positions of all nodes on the side chain of which point is the end node.
# Input : pt is a point object, end node of some side chain
# Output : a list of n+1 node positions. Positions of C-C-C-S for n=3
# Global effect : None
def SideChainNodes(pt, manifolds, nsample):
        result = [pt.pos]
        ind = pt.index
        for i in range(len(manifolds[pt.chain])):
                ind = ind / nsample
                result.append(np.array(manifolds[pt.chain][-(i+1)])[:, ind])
        return result

# NoClash function compares the minimum distance between three side chains and the backbone to the threshold minDist 
# Input : nodesi is the list of array positions of the nodes of side chain i.
# Output : boolean value True if no clash, False otherwise.
# Global effect : None
def NoClash(nodes, Ns, backbone, minDist, minIonDist):
	# comparing all pairwise distances between nodes in the peptoid	
	for node in range(len(nodes)):
		M0 = np.zeros([0,3])
		M0 = np.vstack((M0, nodes[node][-3:]))
		M = np.zeros([0,3])
		for i in [x for x in range(len(nodes)) if x!=node]:
			M = np.vstack((M, nodes[i]))
		M = np.vstack((M, backbone))
		R = scipy.spatial.distance.cdist(M0, M, 'euclidean')
		if (R.min() < minDist):
			return False
	'''
	# comparing the would-be copper location to the backbone
	sum = 0
	for j in range(len(nodes)):
		sum = sum + np.array(nodes[j][-1])
	copper = sum/(float)(len(nodes)) 
	N = np.vstack((backbone, copper))
	S = scipy.spatial.distance.pdist(N, 'euclidean')
	return (nclash and (S.min() > minIonDist))
	'''
	# check that the desirable polygon doesn't cross the backbone
	vertices = []
	for j in range(len(nodes)):
		vertices.append(nodes[j][0])
	tri_bb_index = Choose(len(Ns), 3)
	tri_vertices_index = Choose(len(vertices), 3)
	tri_bb = []
	tri_vertices = []
	for combo in tri_bb_index:
		triangle = []
		for ind in combo:
			triangle.append(Ns[ind])
		tri_bb.append(triangle)
	for combo in tri_vertices_index:
		triangle = []
		for ind in combo:
			triangle.append(list(vertices[ind]))
		tri_vertices.append(triangle)
	for tri1 in tri_bb:
		for tri2 in tri_vertices:
			if Clash(tri1, tri2):
				return False
	return True

# CheckAngle checks that the binding side chains are pointing towards the metal ion to be binded.
# Input : nodesi are lists of array positions of the nodes of side chain i, threshold is the maximum angle difference allowed.
# Output : boolean value, True if the side chains' last bonds are pointing within threshold angle from the S->metalIon angle, False otherwise
# Global effect: None 
def CheckAngle(nodes, threshold):
	sum = 0
	for j in range(len(nodes)):
		sum = sum + np.array(nodes[j][0])
	center = sum/(float)(len(nodes)) 
	for i in range(len(nodes)):
		u = center-nodes[i][0]
		v = nodes[i][0]-nodes[i][1]
		angle = np.arccos(np.dot(u,v)/(LA.norm(u)*LA.norm(v)))
		if angle < threshold:
			continue
		else:
			return False
	return True

# return whether the two triangles are clashing
# Input : tri1, tri2 are lists of three vertices
# Output : true if they clash, false if otherwise
# Global effect: none
def Clash(tri1, tri2):
	#check line segment collision with the triangle tri2
	for i in range(3):
		line = [tri1[i],tri1[(i+1)%3]]
		if LineTriClash(line,tri2):
		      	return True
	return False

# check collision of line with tri
# Input : line is a list of two points, tri is a list of three points
# Output : true if they clash, false if otherwise
# global effect : none
def LineTriClash(line, tri):
	l = [np.array(line[0])-np.array(tri[0]), np.array(line[1])-np.array(tri[0])]
	t = [np.array(tri[1])-np.array(tri[0]), np.array(tri[2])-np.array(tri[0])]
	dl = l[1]-l[0]
	params = np.array([[t[0][0],t[1][0],-dl[0]],[t[0][1],t[1][1],-dl[1]],[t[0][2],t[1][2],-dl[2]]])
	try:
		sol = np.linalg.solve(params,l[0])
	except LinAlgError:
		return False
	if sol[2]<0 or sol[2]>1:
		return False
	else:
		if sol[0]>=0 and sol[1]>=0 and sol[0]+sol[1]<=1:
			return True
	return False

		
