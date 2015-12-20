#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import numpy.random as rd
import scipy.linalg as ln
from HOO import HOO

def POO(f,Nev,K,rhomax,numax):
	Dmax = np.log(K)/np.log(1./rhomax)
	n = 1  #0 dans l'article mais je pense qu'ils etait bourres
	N = 1
	S = [[numax,rhomax]]
	MU = [-np.inf]
	while n < Nev :
		while N >= 0.5*Dmax*np.log(n/np.log(n)) and n<Nev :
			print float(n)/N
			for i in xrange(N) :
				s = [numax,rhomax**(2.*N/(2*i+1))]
				S.append(s)
				T,mu = HOO(f,s[0],s[1],int(float(n)/N))
				MU.append(np.mean(mu))
			n = 2*n
			N = 2*N

		for i,s in enumerate(S) :
			T,mu = HOO(f,s[0],s[1],1)
			MU[i] = MU[i] + np.mean(mu)
		n = n+N
		N = N+1
		sm = S[np.argmax(MU)]
	return sm
