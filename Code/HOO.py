#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import numpy.random as rd
from time import time as ti

def HOO(f,nu,rho,Nev,T=[[0,1]],N=[],mu=[],B=[],U=[],D = {}):
	REW = []
	for n in xrange(1,Nev+1):
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
		H,I = h,i # Seul le dernier noeud n'est pas dans l'arbre.
		a,b = (I-1.)/2.**H,I/2.**H
		X = a + (b-a)*rd.rand()
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
			mu[n_ind] = (1-1./N[n_ind])*mu[n_ind] + Y/N[n_ind]
		for k in xrange(len(T)) :
			h = T[k][0]
			U[k]= mu[k] + np.sqrt(2*np.log(n)/N[k]) + nu*rho**h
		B.append( [np.inf,np.inf] )
		T1 = list(T)
		while T1 != [[0,1]] :
			[h,i] = T1[-1]
			hp,ip = h-1,int((i+1)/2)
			ind = D[str(hp)+str(ip)] #Il faut recuperer l'indice du parent, vu que B[ind(h,i)] contient les b-valeurs des enfants de (h,i)
			B[ind][2*ip-i] = np.min([U[ind],np.max(B[ind])]) # gauche : i = 2*ip, droite: i = 2*ip-1
			T1.remove(T[D[str(h)+str(i)]])
	return T,N,mu,B,U,D,REW # En fait le max est obtenu en x qui maximise mu
