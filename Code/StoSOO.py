#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import numpy.random as rd
import scipy.linalg as ln

def rew(f,node) :
	H,I = node
	alpha,beta = (I-1.)/2.**H,I/2.**H
	X = alpha + (beta-alpha)*rd.rand()
	rew = f(X) # On prend un X dans la boite correspondante et on tire le bras
	return rew

def SOO(f,n,k,hmax,delta):
	T = [[0,1]]	#Initialisation classique
	t = 0 # Nombre d'evaluataions de la fonction effectuees
	MU= [0] # La moyenne des evaluations de la fonction au noeud correspondant dans T (MU[i] = moy des Y pour T[i])
	N = [0] # Nombre d'evaluations faites au noeud T[i]
	break_count = [0]
	b = [np.inf] # B values calculees sur la base de MU + coefficient correctif
	while t < n :
		#print "t=",t
		bmax = -np.inf
		depth = np.max([h for [h,i] in T])
		for h in xrange(min(depth,hmax)+1) :
			B = []
			Lh = [j for j,[p,i] in enumerate(T) if p == h] # tous les (indices des) noeuds de T a profondeur h
			for i in Lh :# on calcule les b-values de tous les noeuds a profondeur h
				if N[i] == 0 :
					b[i] = np.inf # +inf si on y est jamais passe, pour forcer au moins un passage
				else :
					b[i] = MU[i] + np.sqrt(np.log(n*k/delta)/(2*N[i])) # calcul a la UCB
				B.append(b[i])
			im = Lh[np.argmax(B)] # Index of the best node in T at depth h
			if b[im] >= bmax :
				[H,I] = T[im]
				if N[im] < k :# si on a fait moins de k evaluations a ce noeud et que c'est le meilleur, on continue
					r = rew(f,[H,I]) # On prend un X dans la boite correspondante et on tire le bras
					N[im] = N[im] + 1
					MU[im] = (1-1./N[im])*MU[im] + r/N[im]
					t = t+1
				else :
					for s in [[H+1,2*I],[H+1,2*I-1]]:
						if s not in T:
							T.append(s)
							b.append(np.inf)
							MU.append(0)
							N.append(0)
						else :
							MU[im] = -np.inf # Pas ouf de faire ca, mais en theorie c'est la meilleure chose a faire
				bmax = b[im]
	hm = min(depth,hmax) # on va chercher la profondeur maximale atteinte pour laquelle on a calcule toutes les b-values (il se peut que depth = hmax+1 et dans ce cas les dernieres b-values n'ont pas ete calculees)
	MUh = [[MU[j],i] for j,[h,i] in enumerate(T) if h == hm] # Je vais prendre le noeud a profondeur hm qui a la plus grosse MU
	MUh = np.array(MUh)
	im = MUh[np.argmax(MUh[:,0]),1] # et voila je l'ai fait
	return (2*im-1.)/2**(hm+1) # Milieu du meilleur segment
