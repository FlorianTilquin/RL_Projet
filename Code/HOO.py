#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import numpy.random as rd
import scipy.linalg as ln

def HOO(f,Df,nu,rho,P)
	T = [(0,1)]
	B[1,2],B[2,2] = np.inf,np.inf ## Changer les indices, initialiser
	for n = 1,2.
		h,i = 0,1
		P = [(0,1)]
		while (h,i) in T:
			if B[2,2] > B[1,1] :
				h,i = h+1,2*i-1
			if B[2,2] > B[1,1] :
				h,i = h+1,2*i
			else :
				Z = rd.randint(2)
				h,i = h+1,2*i-Z
			P.append((h,i))
		[a,b] = P[h,i]	## Checker que notre partitionning est sous forme de segments.
		X = a + (b-a)*rd.rand()
		Y = f(X)
		T.append((h,i))
		for (h,i) in P:
			T[h,i] += 1
			mu[h,i] = (1-1./T[h,i])*mu[h,i] + Y/T[h,i]
		for (h,i) in T:
			U[h,i]= mu[h,i] +np.sqrt(2*np.log(n)/T[h,i]) + mu*rho**h
		B[h+1,2*i-1],B[h+1,2*i] = np.inf,np.inf
		T1 = T.copy()
		while T1 != [(0,1)]:
			(h,i) = leaf(T1) ##leaf à coder !!
			B[h,i] = np.min([U[h,i],np.max(B[h+1,2*i-1],B[h+1,2*i])])
			T1 = T1.remove((h,i))
	return ??? # En fait il est obtenu où le max ??
