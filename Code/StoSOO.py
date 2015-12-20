#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import numpy.random as rd
import scipy.linalg as ln

def SOO(f,n,k,hmax,delta):
	T = [[0,1]]
	t = 0
	while t < n :
		bmax = -np.inf
		depth = np.max([h for [h,i] in T])
		for h in xrange(min(depth,hmax)+1) :
			B = []
			Lh = [T.index(s) for s in T if S[0] == h]
			for i in Lh :
				b[i] = MU[i] + np.sqrt(np.log(n*k/delta)/(2*N[i]))
				B.append(b[i])
				im = Lh[np.argmax(B)] # Index of the best node in T
			if b[im] >= bmax :
				if N[im] < k :
					[H,I] = T[im]
					a,b = (I-1.)/2.**H,I/2.**H
					X = a + (b-a)*rd.rand()
					rew = f(X)
					t = t+1
				else :
					T.append([H+1,2*I],[H+1,2*I-1])
					MU.extend([0,0])
					bmax = b[im]
	MUh = [MU[i] for i in xrange(len(T)) if T[i][0] == depth]
	muhm = np.min(MUh)
	[h,i] = T[MU.index(muhm))]
	return (2*i-1.)/2**(h+1) # Milieu du meilleur segment
