'''
Copyright (c) 2016, Tian Jiang
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''


def read_residue_carbons(filename):
    CA = []
    CB = []
    CG = []
    firstCB = 1
    firstCG = 1
    with open(filename) as residuefile:
	residuelines = residuefile.readlines()
    for i in range(len(residuelines)):
        splitted = residuelines[i].split()
        if len(splitted) > 8:
            if splitted[2] == "CA":
                CA.append(splitted[6:9])
            if splitted[2] == "CB":
                if firstCB == 1:
                    CB = splitted[6:9]
                    firstCB = 0
                else:
                    assert(CB == splitted[6:9])
            if splitted[2] == "CG":
                if firstCG == 1:
                    CG = splitted[6:9]
                    firstCG = 0
                else:
                    assert(CG == splitted[6:9])
    x = 0.
    y = 0.
    z = 0.
    for pos_str in CA:
        x += float(pos_str[0])
        y += float(pos_str[1])
        z += float(pos_str[2])
    x = x/len(CA)
    y = y/len(CA)
    z = z/len(CA)
    CA = [x,y,z]
    CB = [float(x) for x in CB]
    CG = [float(x) for x in CG]
    return [CA,CB,CG]

def read_backbone_cacb(filename):
    carbons = []
    with open(filename) as bfile:
        blines = bfile.readlines()
    cacb = []
    for i in range(len(blines)):
        if blines[i].startswith('TER'):
            break
        else:
            splitted = blines[i].split()
            if splitted[2] == "CA":
                cacb.append([float(x) for x in splitted[5:8]])
                
            if splitted[2] == "CB":
                cacb.append([float(x) for x in splitted[5:8]])
                carbons.append(cacb)
                cacb = []
    return carbons

if __name__=="__main__":
    #CA = read_residue_carbons("./data/residues/mdm2_leu_stubs_trim.pdb")
    #print CA
    carbons = read_backbone_cacb("./data/backbones/oop_dimer_DDDD.pdb")
    print carbons
