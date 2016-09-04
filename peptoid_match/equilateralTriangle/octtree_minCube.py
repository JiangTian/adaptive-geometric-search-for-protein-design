'''
Copyright (c) 2016, Tian Jiang
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

# octree with minimum cube length
import numpy

class Point:
        def __init__(self, chain, index, position):
            self.chain = chain
	    self.index = index	
            self.pos = position

#note that this order is important. Do not change it
tree_splits = [numpy.array(x, dtype=float) for x in [[ 1, 1, 1], [ 1, 1,-1],
                                                     [ 1,-1, 1], [ 1,-1,-1],
                                                     [-1, 1, 1], [-1, 1,-1],
                                                     [-1,-1, 1], [-1,-1,-1]]]

# points must have a field pos (a 3d array)
class OctTreeBox:
    def __init__(self, center, halfSize, minSize, depth=0):
        self.center = center
        self.halfSize = halfSize
        self.minSize = minSize
        self.isLeaf = True
        self.children = None
        self.points = []

    def SplitLeafToNode(self):
        assert(self.isLeaf)
        qs = self.halfSize / 2
        c = self.center
        new_centers = [numpy.copy(c) + x * qs for x in tree_splits]
        new_boxes = [OctTreeBox(x, qs, self.minSize) for x in new_centers]
        for pt in self.points:
            p = pt.pos
            idx = ((p[0]<c[0]) << 2) + ((p[1]<c[1]) << 1) + (p[2]<c[2])
            new_boxes[idx].AddPoint(pt)
        self.isLeaf = False
        self.points = None
        self.children = new_boxes

    def AddPoint(self, point):
        if self.isLeaf:
            self.points.append(point)
            if self.halfSize >= self.minSize:
                self.SplitLeafToNode()
        else:
            c = self.center
            p = point.pos
            idx = ((p[0]<c[0]) << 2) + ((p[1]<c[1]) << 1) + (p[2]<c[2])
            self.children[idx].AddPoint(point)

    def GetDepth(self):
        #TODO: this could be computed during construction
        if self.isLeaf:
            return 0
        else:
            return 1+max([c.GetDepth() for c in self.children])

    def GetPoints(self):
        if self.isLeaf:
            return self.points
        else:
            #TODO: this is not fast at all...
            return reduce(lambda x, y: x+y, [z.GetPoints() for z in self.children])

    def GetChildren(self):
        # return only non empty children
        assert(not self.isLeaf)
        return [x for x in self.children if x.points or not x.isLeaf]

    def toString(self, indent = 0):
        ind = " " * indent
        if self.isLeaf:
            s  = ind + "Leaf.\n"
            s += ind + "->center: " + str(self.center) + "\n"
            s += ind + "->halfSize: " + str(self.halfSize) + "\n"
            s += ind + "->points:\n"
            for p in self.points:
                s += ind + "  " + str(p.pos) + "\n"
        else:
            s  = ind + "Node.\n"
            s += ind + "->center: " + str(self.center) + "\n"
            s += ind + "->halfSize: " + str(self.halfSize) + "\n"
            s += ind + "->children:\n"
            for c in self.GetChildren():
                s += c.toString(indent+2)
        return s

    def __str__(self):
        return self.toString(0)

def GetOctTree(points, center, halfSize):
    octtree = OctTreeBox(center, halfSize)
    for p in points:
        octtree.AddPoint(p)
    return octtree

if __name__=="__main__":
    # debug
    import numpy.random
    class Pt:
        def __init__(self, p):
            self.pos = p
    N = 100
    center = numpy.array([0.5, 0.5, 0.5])
    halfSize = 0.5
    points = [Pt(x) for x in numpy.random.rand(N, 3)]
    ot = GetOctTree(points, center, halfSize)
    print(ot)
