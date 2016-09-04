'''
Copyright (c) 2016, Tian Jiang
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import os
import shlex
import sys

def writeLine(line,ind):
    [serial, atomName, resName, chainId, resSeq, x, y, z, occupancy, tempFactor] = line
    serial = int(serial)
    resName = 310 # until Doug defines it
    chainId = "A"
    resSeq = ind+1
    x = float(x)
    y = float(y)
    z = float(z)
    occupancy = float(occupancy)
    tempFactor = float(tempFactor)
    return "ATOM  %5d  %-*s %d %s %3d %11.3f %7.3f %7.3f  %.2f  %.2f" %(serial, 3, atomName, resName, chainId, resSeq, x, y, z, occupancy, tempFactor)
# NCCCS
with open("../6mer.pdb") as ba:
    content1 = ba.readlines()
b = [8,13,18,23,28,33]
content2 = [content1[i] for i in b]
for filename in os.listdir(os.getcwd()):
    if filename.endswith(".pdb"):
        with open(filename) as f:
            content = f.readlines()
        bb = [0,4,8]
        content22 = [content[i] for i in bb]
        # returns the list of take-off nitrogen indices
        def check():
            result = []
            for i in range(3):
                sp22 = shlex.split(content22[i])
                c1 = map(float,sp22[6:9])
                for j in range(6):
                    sp2 = shlex.split(content2[j])
                    ca1 = map(float,sp2[6:9])
                    if c1 == ca1:
                        result.append(j)
            return result
        ns = check() # take-off nitrogens
        index = [b[i] for i in ns]
        def insertlines(chain):
            ind = ns[chain]
            if chain == 0:
                a=1
                b=4
            if chain == 1:
                a=5
                b=8
            if chain == 2:
                a=9
                b=12
            for i in range(a,b):
                args = shlex.split(content[i])
                newline = writeLine(args[1:11],ind)
                newfile.write(newline)
                newfile.write('\n')
        
        newfile = open("syn_"+filename,"w")
        for i in range(len(content1)):
            for j in range(3):
                if i == index[j]+1:
                    insertlines(j)
            newfile.write(content1[i])
        newfile.close

