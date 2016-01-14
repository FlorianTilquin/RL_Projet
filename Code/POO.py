#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import numpy.random as rd
import scipy.linalg as ln
from HOO import HOO

def POO(f,Nev,rhomax,numax,noise):
	K = 2.
	Dmax = np.log(K)/np.log(1./rhomax)
	n = 0
	N = 1 # Nombre de HOO lances
	S = [[numax,rhomax]] # Parametres des HOO lances
	MU = [-np.inf]
	H = [[]] #Stock les HOO
	k = 0
	POS_global = []
	REW_global = []
	def condition(N,n):
		if n < 2 :
			return True
		else :
			b = 0.5*Dmax*np.log(n/np.log(n))
			return (N >= b)*(N <= 2*b)
	while n < Nev :
		#print n,N
		#print "cond", condition(n,n)
		while condition(N,n) and n < Nev :
			for i in xrange(N) :
				s = [numax,rhomax**(2.*N/(2*i+1))]
				s = [numax,0.5]
				S.append(s)
				NevH = int(float(n)/N)+1
				dum = None
				#print "Nev",NevH
				#print NevH
				dum = list(HOO(f,NevH,s[0],s[1],noise = noise,T=[],N=[],mu=[],B=[],U=[],D={},Pos=[],REW=[]))
				POS_global.extend(dum[-2])
				REW_global.extend(dum[-1])
				H.append( dum )
				#print "H",len(dum[2]),len(H)
				MU.append(np.max(dum[2])) # On prend le meilleur resultat de cette HOO
			n = n + N * NevH
			N = 2*N
			k = k + NevH
		for i,s in enumerate(S) :
			H[i] = list(HOO(f,1,s[0],s[1],noise,*H[i]))
			POS_global.append(H[i][-2][-1])
			REW_global.append(H[i][-1][-1])
			#Ni = H[i][1]
			MU[i] = np.max(H[i][2]) #MU[i]*(1-1./k) + np.max(H[i][2])/k
		n = n+N
		k = k+1
		sm = np.argmax(MU)
		#print "s,S[s]:",sm,S[sm]
		#print "k", k
	return H[sm],POS_global,REW_global
