#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import numpy.random as rd
import scipy.linalg as ln
from HOO import HOO

def POO(f,Nev,K,rhomax,numax):
	Dmax = np.log(K)/np.log(1./rhomax)
	n = 0 # Nombre d'evaluations de la fonction. Faut pas mettre 0 sinon ca tourne a l'infini et j'ai la flemme de faire un truc qui corrige cette connerie
	N = 1 # Nombre de HOO lances
	S = [[numax,rhomax]] # Parametres des HOO lances
	MU = [-np.inf]
	H = [[]] #Stock les HOO
	k = 0
	def condition(N,n):
		if n < 2 :
			return True
		else :
			b = 0.5*Dmax*np.log(n/np.log(n))
			return (N >= b)*(N <= 2*b)
	while n < Nev :
		print n,N
		while condition(N,n) and n < Nev :
			for i in xrange(N) :
				s = [numax,rhomax**(2.*N/(2*i+1))]
				s = [numax,0.5]
				S.append(s)
				NevH = int(float(n)/N)+1
				dum = list(HOO(f,s[0],s[1],NevH))
				H.append( dum ) # On lance une HOO de n/N evaluations (ce qui en pratique fait pas beaucoup...)
				MU.append(np.max(dum[2])) # On prend le meilleur resultat de cette HOO
			n = n + N * NevH
			N = 2*N
			k = k + NevH
		for i,s in enumerate(S) :
			H[i] = list(HOO(f,s[0],s[1],1,*H[i]))
			#Ni = H[i][1]
			MU[i] = np.max(H[i][2]) #MU[i]*(1-1./k) + np.max(H[i][2])/k
		n = n+N
		k = k+1
		sm = np.argmax(MU)
		#print "s,S[s]:",sm,S[sm]
	return H#list(H[sm])
