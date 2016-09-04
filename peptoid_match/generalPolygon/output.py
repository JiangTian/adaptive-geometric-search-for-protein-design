'''
Copyright (c) 2016, Tian Jiang
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import shlex

# ToDo: draw the full backbone
# Draw the backbone in green, side chains in blue and mark the binding polygon 
def Draw(kpt,nodes,Ns):
	assert(len(kpt) == len(nodes))
	fig = figure()
	ax = Axes3D(fig)
	k = len(kpt)
	x = [[] for i in range(k)]
	y = [[] for i in range(k)]
	z = [[] for i in range(k)]
	px = []
	py = []
	pz = []
	for i in range(k):
		x[i] = np.append(np.array(Ns)[kpt[i].chain,0], np.array(nodes[i])[:,0])
		y[i] = np.append(np.array(Ns)[kpt[i].chain,1], np.array(nodes[i])[:,1])
		z[i] = np.append(np.array(Ns)[kpt[i].chain,2], np.array(nodes[i])[:,2])
		px.append(kpt[i].pos[0])
		py.append(kpt[i].pos[1])
		pz.append(kpt[i].pos[2])
		ax.plot(x[i], y[i], z[i], 'b') #plot the sidechains
	ax.plot(np.append(np.array(Ns)[:,0], np.array(Ns)[0,0]), np.append(np.array(Ns)[:,1], np.array(Ns)[0,1]), np.append(np.array(Ns)[:,2], np.array(Ns)[0,2]),'g') #plot the backbone Ns
	px.append(kpt[0].pos[0])
	py.append(kpt[0].pos[1])
	pz.append(kpt[0].pos[2])
	ax.plot(px, py, pz, 'r', linestyle="--") #plot the desirable polygon
	for j in range(k):
		ax.scatter([np.array(nodes[j])[3,0]], [np.array(nodes[j])[3,1]], [np.array(nodes[j])[3,2]], 'r') #plot the binding nodes
	show()

def writeLine(line):
    [ATOM, serial, atomName, resName, chainId, resSeq, x, y, z, occupancy, tempFactor] = line
    serial = int(serial)
    resName = int(resName)
    resSeq = int(resSeq)
    x = float(x)
    y = float(y)
    z = float(z)
    occupancy = float(occupancy)
    tempFactor = float(tempFactor)
    return "ATOM  %5d  %-*s %d %s %3d %11.3f %7.3f %7.3f  %.2f  %.2f" %(serial, 3, atomName, resName, chainId, resSeq, x, y, z, occupancy, tempFactor)


def writeLink(link):
	positions = [0, 13, 17, 21, 25, 43, 47, 51, 55, 61, 68, 74]
	ret = ""
	for s, p in zip(link, positions):
		while len(ret) < p:
			ret = ret + " "
		ret = ret + s
	ret = ret + "\n"
	return ret

# output a pdb file for this design "nodes" on the backbone
def Output(kpt,nodes,bblines,num):
	chains = []
	filename = '6mer_'
	for p in kpt:
		chains.append(p.chain)
		filename += str(p.chain)
	newfile = open('output/'+filename+'_'+str(num)+'.pdb',"w")
	for i in range(4):
		newfile.write(bblines[i])
	link = shlex.split(bblines[4])
	if 0 in chains:
		link[2] = '310'
	if 5 in chains:
		link[6] = '310'
	newfile.write(writeLink(link))
	a = 0
	atomName = ['CB','CG','SD']
	for i in range(5, len(bblines)-1):
		splited = shlex.split(bblines[i])
		chain_ind = int(splited[5])-1
		if chain_ind in chains:
			a += 1
			splited[3] = '310'
			if a != 5:
				newfile.write(writeLine(splited))
				newfile.write('\n')
			else:   
				a = 0
				newfile.write(writeLine(splited))
				newfile.write('\n')
				cNodes = nodes[chains.index(chain_ind)]
				serial = int(splited[1])
				for j in range(1, len(nodes[0])):
					serial += 1
					line = ['ATOM', serial, atomName[j-1], '310', 'A', chain_ind+1] + list(cNodes[j]) + [1.00, 1.00]
					newfile.write(writeLine(line))
					newfile.write('\n')
		else:
			newfile.write(bblines[i])
	newfile.write(bblines[-1])
	newfile.close
	return 0
