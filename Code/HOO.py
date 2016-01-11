#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import numpy.random as rd
from time import time as ti

def HOO(f,nu,rho,Nev,T=[],N=[],mu=[],B=[],U=[],D = {},Pos =[],REW = []):
	lnev = 2*np.log(Nev)
	#nuro = nu*rho**np.arange(25)
	#Lnev = lnev/np.arange(1,Nev+10)
	for n in xrange(1,Nev+1):
		#print "n=",n
		#if n %50 == 0:
		#	print n
		P = [[0,1]]
		h,i = 0,1
		while [h,i] in T:
			j = D[str(h)+str(i)]
			[b1,b2] = B[j]
			[h,i] = T[j]
			if b2 > b1 :
				h,i = h+1,2*i-1
			if b2 < b1 :
				h,i = h+1,2*i
			if b2 == b1 :
				Z = rd.randint(2)
				h = h+1
				i = 2*i - Z
			P.append([h,i])
		H,I = h,i
		a,b = (I-1.)/2.**H,I/2.**H
		X = a + (b-a)*rd.rand()
		Pos.append(X)
		Y = f(X)
		REW.append(Y)
		T.append([H,I])
		D[str(H)+str(I)] = len(T)-1
		N.append(1)
		mu.append(Y)
		U.append(mu)
		for [k,l] in P :
			n_ind = D[str(k)+str(l)]
			N[n_ind] += 1
			mu[n_ind] *= (1-1./N[n_ind])
			mu[n_ind] += Y/N[n_ind]
			U[n_ind] = mu[n_ind] + lnev/N[n_ind] + nu*rho**k
		B.append( [np.inf,np.inf] )
		P1 = list(P)
		while P1 != [[0,1]] :
			[h,i] = P1[-1]
			hp,ip = h-1,(i+1)/2
			ind = D[str(hp)+str(ip)]
			B[ind][2*ip-i] = np.min([U[ind],np.max(B[ind])]) # gauche : i = 2*ip, droite: i = 2*ip-1
			P1.pop()
	return T,N,mu,B,U,D,Pos,REW
