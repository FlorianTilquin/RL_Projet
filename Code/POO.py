#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import numpy.random as rd
import scipy.linalg as ln
from HOO import HOO

def POO(f,Nev,K,rhomax,numax):
	Dmax = np.log(K)/np.log(1./rhomax)
	n = 1 # Nombre d'evaluations de la fonction. Faut pas mettre 0 sinon ca tourne a l'infini et j'ai la flemme de faire un truc qui corrige cette connerie
	N = 1 # Nombre de HOO lances
	S = [[numax,rhomax]] # Parametres des HOO lances
	MU = [-np.inf]
	H = [[]]
	k = 0
	def borne(n) :
		if n < 2 :
			return 0
		else :
			return 0.5*Dmax*np.log(n/np.log(n))
	while n < Nev :
		while N >= borne(n)  and n < Nev : # j'ai un doute sur la capacite de python a se rendre compte de la valeur de n en permance, je ferai des tests c'est interessant. au passage note la pertinence du truc si n=0 ou 1...
			for i in xrange(N) : # on va lancer N HOO sur une grille exponentielle
				s = [numax,rhomax**(2.*N/(2*i+1))]# ce que j'appelle "grille exponentielle"
				S.append(s)
				NevH = int(float(n)/N)+1
				H.append( list(HOO(f,s[0],s[1],NevH)) ) # On lance une HOO de n/N evaluations (ce qui en pratique fait pas beaucoup...)
				MU.append(np.max(H[-1][2])) # la je sais pas trop quoi faire d'aautre que la moyenne... Le meilleur resultat sinon ?alllezz j'essaie
			n = n + N * NevH
			N = 2*N
			k = k + NevH
		for i,s in enumerate(S) :
			H[i] = list(HOO(f,s[0],s[1],1,*H[i]))
			Ni = H[i][1]
			MU[i] = MU[i]*(1-1./k) + np.max(H[i][2])/k
		n = n+N
		k = k+1
		sm = np.argmax(MU)
	return H[sm]
