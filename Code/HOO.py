#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import numpy.random as rd
import scipy.linalg as ln

def HOO(f,nu,rho,Nev):
	#T = [[0,1]]
	#N = [1]
	#mu = [-0.5]
	#B = [[np.inf,np.inf]]
	#U = [0.]
	T,N,mu,B,U = [],[],[],[],[]
	for n in range(1,Nev+1):
		print( str(n)+" evaluations de la fonction effectuees" )
		P = [[0,1]]
		h,i = 0,1
		while [h,i] in T:
			j = T.index([h,i])
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
		H,I = h,i # seul le dernier noeud n'est pas dans l'arbre.
		a,b = (I-1.)/2.**H,I/2.**H
		X = a + (b-a)*rd.rand()
		Y = f(X)
		T.append([H,I])
		N.append(1)
		mu.append(Y)
		U.append(mu)
		for [k,l] in P :
			n_ind = T.index([k,l])
			N[n_ind] += 1
			mu[n_ind] = (1-1./N[n_ind])*mu[n_ind] + Y/N[n_ind]
		for k in range(len(T)) :
			h = T[k][0]
			U[k]= mu[k] + np.sqrt(2*np.log(n)/N[k]) + nu*rho**h
		B.append( [np.inf,np.inf] )
		T1 = list(T)
		while T1 != [[0,1]] :
			[h,i] = T1[-1]
			hp,ip = h-1,int((i+1)/2)
			#print T
			#print P
			#print h,i
			ind = T.index([hp,ip]) #Il faut recuperer l'indice du parent, vu que B[ind(h,i)] contient les b-valeurs des enfants de (h,i)
			gd = 2*ip-i # C'est la branche de gauche ou droite ?gauche : i =2*ip, droite: i = 2*ip-1
			B[ind][gd] = np.min([U[ind],np.max(B[ind])])
			T1.remove([h,i])
	return T,N,mu # En fait il est obtenu ou le max ??
