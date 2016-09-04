'''
Copyright (c) 2016, Tian Jiang
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import numpy as np

# return a list of (n choose k) lists
def Choose(n,k):
	if n==k:
		return [[x for x in range(n)]]
	elif k==0:
		return [[]]
	else:
		result = []
		lis = Choose(n-1,k-1)
		for i in range(len(lis)):
			result.append([0]+list(np.array(lis[i])+1))
		lit = Choose(n-1,k)
		for i in range(len(lit)):
			result.append(list(np.array(lit[i])+1))
		return result

# returns all arrangements of the elements of s(set) in a list
def Arrange(S):
	if len(S) == 1:
		return [[S.pop()]]
	else:
		result = []
		for s in S:
			S0 = set(S)
			S0.remove(s)
			res0 = Arrange(S0)
			for l in res0:
				result.append([s]+l)
		return result

if __name__=="__main__":
	n=5
	k=2
	p = Arrange(set([1,2,3,4]))
	print p
