#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import numpy.random as rd
import scipy.linalg as ln

def SOO(f,n,k,hmax,delta):
	T = [[0,1]]
	t = 0
	MU= [0]
	N = [0]
	b = [np.inf]
	print np.sqrt(np.log(n*k/delta)/(2*1.))
	while t < n :
		bmax = -np.inf
		depth = np.max([h for [h,i] in T])
		for h in xrange(min(depth,hmax)+1) :
			B = []
			Lh = [ j for j,[p,i] in enumerate(T) if p == h] # tous les noeuds de T a prof h
			for i in Lh :
				if N[i] == 0 :
					b[i] = np.inf
				else :
					b[i] = MU[i] + np.sqrt(np.log(n*k/delta)/(2*N[i]))
				B.append(b[i])
				im = Lh[np.argmax(B)] # Index of the best node in T at depth h
			if b[im] >= bmax :
				[H,I] = T[im]
				if N[im] < k :
					alpha,beta = (I-1.)/2.**H,I/2.**H
					X = alpha + (beta-alpha)*rd.rand()
					rew = f(X)
					if N[im] == 0 :
						MU[im] = rew
						t = t+1
						N[im] = N[im] + 1
					else :
						MU[im] = (1-1./N[im])*MU[im] + rew/N[im]
						t = t+1
						N[im] = N[im] + 1
				else :
					for s in [[H+1,2*I],[H+1,2*I-1]]:
						if s not in T:
							T.append(s)
							b.append(np.inf)
							MU.append(0)
							N.append(0)
				bmax = b[im]
		#print b
	hm = min(depth,hmax)
	MUh = [[MU[j],i] for j,[h,i] in enumerate(T) if h == hm]
	print MUh
	MUh = np.array(MUh)
	im = MUh[np.argmax(MUh[:,0]),1]
	print hm,im
	return (2*im-1.)/2**(hm+1) # Milieu du meilleur segment
