'''
Copyright (c) 2016, Tian Jiang
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import shlex
import numpy as np
import re
import os

def Draw(CAs, CBs, residue, CAs_res):
	fig = figure()
	ax = Axes3D(fig)
	for i in range(3):
		x = [CAs[i][0],CBs[i][0]]
		y = [CAs[i][1],CBs[i][1]]
		z = [CAs[i][2],CBs[i][2]]
		ax.plot(x,y,z,"b")
	xcb = []
	ycb = []
	zcb = []
	for i in range(3):
		xcb.append(CBs[i][0])
		ycb.append(CBs[i][1])
		zcb.append(CBs[i][2])
	xcb.append(CBs[0][0])
	ycb.append(CBs[0][1])
	zcb.append(CBs[0][2])
	ax.plot(xcb,ycb,zcb,"r", linestyle="--")
	for i in range(3):
		x = [CAs_res[i][0],residue[i][1][0],residue[i][2][0]]
		y = [CAs_res[i][1],residue[i][1][1],residue[i][2][1]]
		z = [CAs_res[i][2],residue[i][1][2],residue[i][2][2]]
		ax.plot(x,y,z,"g")
	show()
		

def writeLine(line):
    [ATOM, serial, atomName, resName, chainId,  x, y, z, occupancy, tempFactor] = line
    serial = int(serial)
    chainId = int(chainId)
    x = float(x)
    y = float(y)
    z = float(z)
    occupancy = float(occupancy)
    tempFactor = float(tempFactor)
    return "ATOM  %5d  %-*s %s %5d %11.3f %7.3f %7.3f  %.2f  %.2f" %(serial, 3, atomName, resName, chainId, x, y, z, occupancy, tempFactor)

def writeLine2(line):
    [ATOM, serial, atomName, resName, resID, chainId,  x, y, z, occupancy, tempFactor] = line
    serial = int(serial)
    chainId = int(chainId)
    x = float(x)
    y = float(y)
    z = float(z)
    occupancy = float(occupancy)
    tempFactor = float(tempFactor)
    return "ATOM  %5d  %-*s %s %s %3d %11.3f %7.3f %7.3f  %.2f  %.2f" %(serial, 3, atomName, resName, resID, chainId, x, y, z, occupancy, tempFactor)

# output a pdb file for this design "nodes" on the backbone
def Output(bbfile, ind, M,t1,t2,M_res,residue):
	with open("./data/backbones/"+bbfile) as bfile:
		blines = bfile.readlines()
	ind_str = "".join([str(x) for x in ind])
	newfile = open('output/'+bbfile+'_'+ind_str+'.pdb',"w")
	copy = 0
	for line in blines:
		if line.startswith("TER"):
			newfile.write(line)
			copy = 1
		elif copy == 1:
			newfile.write(line)
		else:
			splitted = line.split()
			if re.search("H",splitted[2]):
				continue
			else:
				v = [float(x) for x in splitted[5:8]]
				w = np.dot(M, np.array(v)-t1)+t2
				w = np.squeeze(np.asarray(w))
				for i in range(3):
					splitted[i+5] = str(w[i])
				newfile.write(writeLine(splitted))
				newfile.write('\n')
	newfile.close

	# use M_res=[M1,M2,M3] to transform the first 5 atoms in each residue and use writeLine2 
	# TODO : rotate the C, O, etc after CA on residues
	# TODO : in residues_uni, I'm taking the first model of each residue but ideally should take the average of all models.
	reslist = sorted(os.listdir("./data/residues_uni"))
	for j in range(len(reslist)):
		resfile = reslist[j]
		t = residue[j][1] # the order of residue must match with M_res and sorted os.listdir("./data/residues")
		with open("./data/residues_uni/"+resfile) as rfile:
			rlines = rfile.readlines()
      		nresfile = open('output/'+resfile,'w') 
		for i in range(len(rlines)):
			splitted = rlines[i].split()
			if i < 5:
				v = [float(x) for x in splitted[6:9]]
				w = np.dot(M_res[j], np.array(v)-t)+t
				w = np.squeeze(np.asarray(w))
				for i in range(3):
					splitted[i+6] = str(w[i])
				nresfile.write(writeLine2(splitted))
				nresfile.write('\n')
			elif len(splitted)>2 and re.search("H", splitted[2]):
				continue
			else:
				nresfile.write(rlines[i])
		nresfile.close
	return 0
