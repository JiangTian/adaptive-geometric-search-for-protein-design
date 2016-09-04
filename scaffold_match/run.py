'''
Copyright (c) 2016, Tian Jiang
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

# scaffold matcher for oop backbones with hot residues binding with mdm2. Default setting gives the best two matching results. 

import numpy as np
from numpy import linalg as LA
import os
import time
from choose import *
from get_data import *
from output import *

starttime = time.time()
########################## Params #######################################
kgon = 3 # for this problem there are three residues
# set error bounds:
# err = 0.1
# angle_err = 0.1 # for the best matching result
#err = 0.3 # maximum tolerance for positions match
#angle_err = 0.4 # maximum tolerance for bond angles match
err = 100000.
angle_err = 100000.
nsample = 100
print 'err = ', err
print 'angle_err = ', angle_err

############################# Read in data ###########################\
residue_files = sorted(os.listdir("./data/residues/"))
assert(len(residue_files) == kgon)
residue = []
for i in range(kgon):
	residue.append(read_residue_carbons("./data/residues/"+residue_files[i])) # reads [ca,cb,cg] for each residue
CB_res = []
for i in range(kgon):
	CB_res.append(residue[i][1])
backbone_c = []
backbone_list = []
for backbonefile in os.listdir("./data/backbones/"):
	if backbonefile.startswith("oop_"):
		backbone_list.append(backbonefile)
		backbone_c.append(read_backbone_cacb("./data/backbones/"+backbonefile)) # reads [ca,cb] for each backbone

##################### auxiliary functions #####################
# computes the side lengths given a triangle
def lengths(tri):
	length = [LA.norm(np.array(tri[0])-np.array(tri[1])), LA.norm(np.array(tri[1])-np.array(tri[2])), LA.norm(np.array(tri[2])-np.array(tri[0]))]
	return length

# check if two triangles have the same shape
def match(tri1,tri2):
	ind = []
	for l1 in tri1:
		for i in range(len(tri2)):
			l2 = tri2[i]
			if l1 <= l2+err and l1 >= l2-err:
				ind.append(i)
	return ind

# returns the common vertex index shared by triangle sides a and b.
def vertice(a,b):
	if set([a,b]) == set([1,2]):
		return 2
	if set([a,b]) == set([1,0]):
		return 1
	if set([a,b]) == set([0,2]):
		return 0

# check if any of the 2-combos in CBs, the list of CB atoms from the backbone, matches with tri, the triangle of CBs on the 3 residues
# returns a list of indices that match with CB_res, in the order of matching vertices with CB_res vertices
def triangle_match(CBs,tri):
	tri1 = lengths(tri)
	inds = Choose(len(CBs),3)
	result = []
	for i in range(len(inds)):
		triangle = [CBs[x] for x in inds[i]]
		tri2 = lengths(triangle)
		tri2_ind = match(tri1,tri2) 
		if set(tri2_ind) == set([0,1,2]):
			vertices_ind = [None,None,None]
			vertices_ind[0] = vertice(tri2_ind[2],tri2_ind[0])
			vertices_ind[1] = vertice(tri2_ind[0],tri2_ind[1])
			vertices_ind[2] = vertice(tri2_ind[1],tri2_ind[2])
			result.append([inds[i][x] for x in vertices_ind])      
	return result

# returns 3D rotational matrix of angle theta around u
def rotate(u, theta):
	M = np.matrix([[np.cos(theta)+u[0]**2*(1-np.cos(theta)), u[0]*u[1]*(1-np.cos(theta))-u[2]*np.sin(theta), u[0]*u[2]*(1-np.cos(theta)) + u[1]*np.sin(theta)],
	[u[0]*u[1]*(1-np.cos(theta)) + u[2]*np.sin(theta), np.cos(theta) + u[1]**2*(1-np.cos(theta)), u[1]*u[2]*(1-np.cos(theta))  - u[0]*np.sin(theta)],
	[u[2]*u[0]*(1-np.cos(theta)) - u[1]*np.sin(theta), u[2]*u[1]*(1-np.cos(theta)) + u[0]*np.sin(theta), np.cos(theta)+u[2]**2*(1-np.cos(theta))]])
	return M

# returns a rotation matrix that rotates u onto v, both are np.arrays
def rot_uv(u,v):
	assert(abs(LA.norm(u) - LA.norm(v)) <= err)
	vec = np.cross(u,v)
	uvec = vec/LA.norm(vec)
	theta = np.arccos(np.dot(u,v)/(LA.norm(u)*LA.norm(v)))
	M = rotate(uvec, theta)
	assert(LA.norm(np.squeeze(np.asarray(M*np.matrix(u).transpose()))-v) <= err)
	return M

# check if residue CA's can match with backbone CA's given matched CBs
def CA_match(CB_bb, CA_bb, CB_res, residue):
	t1 = np.array(CB_bb[0])
	t2 = np.array(CB_res[0]) 
	target = [np.array(CB_res[i]) - np.array(CB_res[0]) for i in range(3)]
	tomove = [np.array(CB_bb[i]) - np.array(CB_bb[0]) for i in range(3)]
	M1 = rot_uv(tomove[1],target[1])
	tomove_int = M1 * np.matrix(tomove[2]).transpose()
	tomove_array = np.squeeze(np.asarray(tomove_int))
	n1 = np.cross(target[1], tomove_array)
	n2 = np.cross(target[1], target[2])
	M2 = rot_uv(n1/LA.norm(n1),n2/LA.norm(n2))
	M = M2 * M1
	new_CA_bb = [np.squeeze(np.asarray(M*np.matrix(np.array(a)-t1).transpose())) + t2 for a in CA_bb]
	new_CB_bb = [np.squeeze(np.asarray(M*np.matrix(np.array(a)-t1).transpose())) + t2 for a in CB_bb]
	anglerr = 0
	for i in range(3):
		cc1 = np.array(residue[i][0]) - np.array(CB_res[i])
		cc2 = np.array(residue[i][2]) - np.array(CB_res[i])
		c1 = np.array(new_CA_bb[i]) - np.array(CB_res[i])
		theta1 = np.arccos(np.dot(cc1,cc2)/(LA.norm(cc1)*LA.norm(cc2)))
		theta2 = np.arccos(np.dot(c1,cc2)/(LA.norm(c1)*LA.norm(cc2)))
		anglediff = theta1 - theta2
		anglerr += anglediff**2
		if abs(theta1-theta2) > angle_err:
			print "Angle mismatch by ", theta1-theta2
			return None, None, None,None,None
	anglerr = (anglerr/3.)**.5
	return M, t1,t2,new_CA_bb, new_CB_bb, anglerr

# gives the closest CA position on residue while maintaining the bond angles
def rotate_res(CAs,residue,nsample):
	CA_res = [0,0,0]
	trans = [0,0,0]
	for i in range(3):
		res = residue[i]
		u = np.array(res[2])-np.array(res[1])
		u = u/LA.norm(u)
		Phi = np.arange(0.0, 2*np.pi, 2*np.pi/nsample)
		minl = 100.
		for phi in Phi:
			M = rotate(u,phi)
			ca = np.dot(M,np.array(res[0])-np.array(res[1])) + np.array(res[1])
			ca = np.squeeze(np.asarray(ca))
			if LA.norm(ca - CAs[i]) < minl:
				minl = LA.norm(ca-CAs[i])
				CA_res[i] = ca
				trans[i] = M
	return CA_res,trans

###################### auxiliary functions for statistic graphs ##################

def match_cal(tri1, tri2):
	minind = [0,1,2]
	minerr = 10000000.
	for tri in  Arrange(tri2):
		norm_temp = LA.norm(np.array(tri1) - np.array(tri))
		if norm_temp < minerr:
			minerr = norm_temp
			minind = [tri2.index(tri[0]), tri2.index(tri[1]), tri2.index(tri[2])]
	return minind, minerr/np.sqrt(3)

def triangle_cal(CBs,tri):
	tri1 = lengths(tri)
	inds = Choose(len(CBs),3)
	result = []
	errs = []
	vertices_ind = [None,None,None]
	for i in range(len(inds)):
		triangle = [CBs[x] for x in inds[i]]
		tri2 = lengths(triangle)
		tri2_ind, mserr = match_cal(tri1,tri2) 
		assert(set(tri2_ind) == set([0,1,2]))
		vertices_ind[0] = vertice(tri2_ind[2],tri2_ind[0])
	      	vertices_ind[1] = vertice(tri2_ind[0],tri2_ind[1])
	       	vertices_ind[2] = vertice(tri2_ind[1],tri2_ind[2])
	       	result.append([inds[i][x] for x in vertices_ind])
		errs.append(mserr)
	return result, errs

			
##################### main loop : search #################
mserrs_list = []
anglerrs_list = []
for backbone in backbone_c:
	#print "---------------"
	CB_backb = []
	for cc in backbone:
		CB_backb.append(cc[1])
	#print "backbone CBs ", CB_backb
	#print "residue_CBs", CB_res
	##inds = triangle_match(CB_backb,CB_res) # this function can be replaced by adaptive search for longer residues
	inds, mserrs = triangle_cal(CB_backb, CB_res) 
	if len(inds) == 0:
		continue # no CB combos from this backbone matches, so continue
	# otherwise check if CAs can match
        for i in range(len(inds)):
	        ind = inds[i]
	        mserr = mserrs[i]
		CB_bb = [CB_backb[x] for x in ind] # matching CBs in the matched order
		CA_bb = [backbone[x][0] for x in ind]# corresponding CAs in order
		trans_bb, t1,t2,CAs, CBs, anglemserr = CA_match(CB_bb, CA_bb, CB_res, residue)#rotate CB_bb onto CB_res and check whether CAs' can match.
		#print mserr
		#print anglemserr
		mserrs_list.append(mserr)
		anglerrs_list.append(anglemserr)
		if trans_bb != None:
			bb_index = backbone_c.index(backbone)
			bbs = [backbone[x] for x in ind]
			CAs_res,trans_res = rotate_res(CAs,residue,nsample)
 			#Draw(CAs,CBs,residue, CAs_res)
			#Output(backbone_list[bb_index],ind,trans_bb,t1,t2,trans_res,residue)

newerrs = sorted(mserrs_list)
for nerr in newerrs:
	print nerr
print "----------"
newangles = sorted(range(len(mserrs_list)), key=lambda k: mserrs_list[k])
for angle_ind in newangles:
	print anglerrs_list[angle_ind]
print "total time used is %.2f seconds" %(time.time() - starttime)
