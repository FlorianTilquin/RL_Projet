#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import numpy.random as rd
import scipy.linalg as ln

def HOO(f,Df,nu,rho,Nev)
	T = [[0,1]]
	N = [1]
	mu = [0]
	B = [[np.inf,np.inf]]
	U = [0]
	for n in xrange(Nev)
		h,i = 0,1
		P = [[0,1]]
		for j in xrange(len(T)):
			[b1,b2] = B[j]
			[h,i] = T[j]
			if b2 > b1 :
				h,i = h+1,2*i-1
			if b2 > b1 :
				h,i = h+1,2*i
			else :
				Z = rd.randint(2)
				h,i = h+1,2*i-Z
			P.append([h,i])
		H,I = h,i
		a,b = (I-1)/2**(H-1),i/2**(H-1)
		X = a + (b-a)*rd.rand()
		Y = f(X)
		if [H,I] not in T:
			T.append([H,I])
			N.append(1)
			mu.append(Y)
		for [k,l] in P:
			n_ind = T.index([k,l])
			N[n_ind] += 1
			mu[n_ind] = (1-1./N[n_ind])*mu[n_ind] + Y/N[n_ind]
		for k in xrange(len(T)):
			h = T[k][0]
			U[k]= mu[k] +np.sqrt(2*np.log(n)/N[k]) + nu*rho**h
		B.append( [np.inf,np.inf] )
		T1 = T.copy()
		while T1 != [[0,1]]:
			[h,i] = T1[-1] ## est-ce qu'on peut vraiment partir de la fin et remonter ?
			ind = T.index([h,i]
			B[ind] = np.min([U[ind],np.max(B[ind])])
			T1 = T1.remove([h,i])
	return ??? # En fait il est obtenu où le max ??
