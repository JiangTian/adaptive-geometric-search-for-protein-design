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
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from octtree_minCube import *
from output import *

#TODO: using graphs notation to simplify the loops in SearchCycles.
def Graph(pairs):
	graph = {}
	for (i,j) in pairs:
		if i not in graph:
			graph[i] = []
		graph[i] += [j]
	return graph

def SearchGraph(nPairs):	
	result = []
	g12 = Graph(nPairs[0])
	g23 = Graph(nPairs[1])
	g31 = Graph(nPairs[2])
	set1 = set([p for (p,q) in nPairs[0]])
	for c1 in set1:
		for c2 in g12[c1]:
			if c2 in g23:
				for c3 in g23[c2]:
					if (c3 in g31) and (c1 in g31[c3]):
						assert(c1.isLeaf and c2.isLeaf and c3.isLeaf)
						result.append([c1,c2,c3])
	return result

# Function SearchCyles loops through all possibilities and produces n cycles
# Input : Npairs is a list of 3(for triangles) lists of 2 tuples.
# Output : a list of triangles
# Global effect : None
def SearchCyclesApprox(nPairs, manifolds, nsample, minDist, minIonDist, Ns, backbonelines, ncycles=None, draw = None, output=None):
	result = []
	Cs = SearchGraph(nPairs)
	for [c1,c2,c3] in Cs:
		for p1 in c1.GetPoints():
			for p2 in c2.GetPoints():
				for p3 in c3.GetPoints():
					nodes1 = SideChainNodes(p1, manifolds, nsample)					
					nodes2 = SideChainNodes(p2, manifolds, nsample)
					nodes3 = SideChainNodes(p3, manifolds, nsample)
					if NoClash(nodes1, nodes2, nodes3, Ns, minDist, minIonDist) and CheckAngle(nodes1, nodes2, nodes3, np.pi/3):		
						result.append([p1,p2,p3])
						if draw:
							Draw(p1,p2,p3,nodes1,nodes2,nodes3,Ns)
						if output:
							Output([p1,p2,p3],[nodes1,nodes2,nodes3],backbonelines,len(result))
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
        #ind2 = pt.index / nsample
	#ind1 = ind2 / nsample
	#return [np.array(manifolds[pt.chain][0])[:, 0], np.array(manifolds[pt.chain][1])[:, ind1], np.array(manifolds[pt.chain][2])[:, ind2], pt.pos]

# NoClash function compares the minimum distance between three side chains and the backbone Ns to the threshold minDist 
# Input : nodesi is the list of array positions of the nodes of side chain i.
# Output : boolean value True if no clash, False otherwise.
# Global effect : None
def NoClash(nodes1, nodes2, nodes3, Ns, minDist, minIonDist):
	truth = False
	# comparing all pairwise distances between nodes in the peptoid	
	M = np.vstack((nodes1, nodes2, nodes3, Ns))
	R = scipy.spatial.distance.pdist(M, 'euclidean')
	nclash = (R.min() > minDist)
	# comparing the would-be copper location to the Ns 
	copper = (nodes1[-1]+nodes2[-1]+nodes3[-1])/3.0 
	N = np.vstack((Ns, copper))
	S = scipy.spatial.distance.pdist(N, 'euclidean')
	return (nclash and (S.min() > minIonDist))
	

# CheckAngle checks that the binding side chains are pointing towards the metal ion to be binded.
# Input : nodesi are lists of array positions of the nodes of side chain i, threshold is the maximum angle difference allowed.
# Output : boolean value, True if the side chains' last bonds are pointing within threshold angle from the S->metalIon angle, False otherwise
# Global effect: None 
def CheckAngle(nodes1, nodes2, nodes3, threshold):
	center = (nodes1[0]+nodes2[0]+nodes3[0])/3.0
	u = center-nodes1[0]
	v = nodes1[0]-nodes1[1]
	angle1 = np.arccos(np.dot(u,v)/(LA.norm(u)*LA.norm(v)))
	u = center-nodes2[0]
	v = nodes2[0]-nodes2[1]
	angle2 = np.arccos(np.dot(u,v)/(LA.norm(u)*LA.norm(v)))
	u = center-nodes3[0]
	v = nodes3[0]-nodes3[1]
	angle3 = np.arccos(np.dot(u,v)/(LA.norm(u)*LA.norm(v)))
	if max(angle1,angle2,angle3) < threshold:
		res = True
	else:
		res = False
	return res

# Draw the backbone in green, side chains in blue and mark the binding polygon
def Draw(p1,p2,p3,nodes1,nodes2,nodes3,Ns):
	fig = figure()
	ax = Axes3D(fig)
	ax.plot(np.append(np.array(Ns)[p1.chain,0], np.array(nodes1)[:,0]), np.append(np.array(Ns)[p1.chain,1], np.array(nodes1)[:,1]), np.append(np.array(Ns)[p1.chain,2], np.array(nodes1)[:,2]), 'b')
	ax.plot(np.append(np.array(Ns)[p2.chain,0], np.array(nodes2)[:,0]), np.append(np.array(Ns)[p2.chain,1], np.array(nodes2)[:,1]), np.append(np.array(Ns)[p2.chain,2], np.array(nodes2)[:,2]), 'b')
	ax.plot(np.append(np.array(Ns)[p3.chain,0], np.array(nodes3)[:,0]), np.append(np.array(Ns)[p3.chain,1], np.array(nodes3)[:,1]), np.append(np.array(Ns)[p3.chain,2], np.array(nodes3)[:,2]), 'b')
	ax.plot(np.append(np.array(Ns)[:,0], np.array(Ns)[0,0]), np.append(np.array(Ns)[:,1], np.array(Ns)[0,1]), np.append(np.array(Ns)[:,2], np.array(Ns)[0,2]),'g')
	ax.plot([np.array(nodes1)[3,0],np.array(nodes2)[3,0],np.array(nodes3)[3,0],np.array(nodes1)[3,0]],
		[np.array(nodes1)[3,1],np.array(nodes2)[3,1],np.array(nodes3)[3,1],np.array(nodes1)[3,1]],
		[np.array(nodes1)[3,2],np.array(nodes2)[3,2],np.array(nodes3)[3,2],np.array(nodes1)[3,2]], 'r', linestyle="--")
	ax.scatter([np.array(nodes1)[3,0]], [np.array(nodes1)[3,1]], [np.array(nodes1)[3,2]], 'r')
	ax.scatter([np.array(nodes2)[3,0]], [np.array(nodes2)[3,1]], [np.array(nodes2)[3,2]], 'r')
	ax.scatter([np.array(nodes3)[3,0]], [np.array(nodes3)[3,1]], [np.array(nodes3)[3,2]], 'r')
	show()
