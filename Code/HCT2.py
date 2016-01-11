#! /usr/bin/env python

import numpy as np
import numpy.random as rd

def HCT(f,Nev,nu,rho,delta):
	c = 2*np.sqrt(1./(1-rho))
	c1 = (rho/(3*nu))**(0.125)
	par = Nev,c,c1,delta,nu,rho
	t = 1
	T = [[0,1],[1,1],[1,2]]
	D = {}
	D['01'],D['11'],D['12'] = 0,1,2
	H = 1
	U = [0, np.inf, np.inf]
	B = [0, np.inf, np.inf]
	N = [0,0,0]
	Mu = [0.,0.,0.]
	POS = []
	REW = []
	while t < Nev+1 :
		#print t
		if np.log2(t) == np.floor(np.log2(t)):
			print "alignement des 12 lunes"
			for i in xrange(len(T)) :
				if N[i] > 0:
					U[i] = Mu[i] + nu*rho**T[i][0] + c*np.sqrt( np.log(1./min(c1*delta/t,1.))/N[i])
			for h in xrange(H,-1,-1):
				Lh = [j for j,[p,i] in enumerate(T) if p == h]
				for i in Lh :
					hp,ip = T[i]
					if [hp+1,2*ip] not in T :
						B[i] = U[i]
					else :
						ichild = D[str(hp+1)+str(2*ip-1)]
						B[i] = min(U[i],max(B[ichild],B[ichild+1]))
		[h,i],P = OptTraverse(T,B,D,N,t,par)
		alpha,beta = (i-1.)/2.**h,i/2.**h
		x = rd.uniform(alpha,beta)
		POS.append(x)
		r = f(x)
		REW.append(r)
		t = t + 1
		tp = 2**(np.floor(np.log2(t))+1)
		for [h,i] in P:
			idx = D[str(h)+str(i)]
			N[idx] += 1
			Mu[idx] = (1.-1./N[idx])*Mu[idx] +r/N[idx]
			U[idx] = Mu[idx] + nu*rho**T[idx][0] + c*np.sqrt( np.log(1./min(c1*delta/tp,1.))/N[idx])
		B = UpdateB(T,P,[h,i],D,U,B)
		tau = c**2*np.log(1/min(c1*delta/tp,1.))/(nu**2*rho**(2*h))
		#tau = c**2/(nu**2*rho**(2*h))
		if N[idx] >= tau and [h+1,2*i] not in T :
			print h
			D[str(h+1)+str(2*i-1)] = len(T)
			D[str(h+1)+str(2*i)] = len(T)+1
			T.extend([[h+1,2*i-1],[h+1,2*i]])
			U.extend([np.inf,np.inf])
			N.extend([0,0])
			Mu.extend([0,0])
			B.extend([np.inf,np.inf])
			if h+1 > H:
				H += 1
	return	T,N,Mu,U,B,D,POS,REW

def OptTraverse(T,B,D,N,t,par):
	tp = 2**(np.floor(np.log2(t))+1)
	Nev,c,c1,delta,nu,rho = par
	h,i = 0,1
	P = [[0,1]]
	idx = 0
	tau = 1. #le vrai tau ?N[0] ?
	while [h+1,2*i] in T and N[idx] >= tau:
		#print "zbra"
		idx = D[str(h+1)+str(2*i-1)]
		if B[idx]>B[idx+1]:
			h,i = h+1,2*i-1
		elif B[idx]<B[idx+1]:
			h,i = h+1,2*i
		else :
			Z = rd.randint(2)
			h,i = h+1,2*i-Z
		P.append([h,i])
		tau = c**2*np.log(1/min(c1*delta/tp,1.))/(nu**2*rho**(2*h))
		#tau = c**2/(nu**2*rho**(2*h))
	return [h,i],P

def UpdateB(T,P,node,D,U,B):
	[h,i] = node
	idx = D[str(h)+str(i)]
	if [h+1,2*i] not in T :
		B[idx] = U[idx]
	else :
		ichild = D[str(h+1)+str(2*i-1)]
		B[idx] = min(U[idx],max(B[ichild],B[ichild+1]))
	P.pop()
	for [h,i] in P[-2::-1]:
		idx = D[str(h)+str(i)]
		ichild = D[str(h+1)+str(2*i-1)]
		B[idx] = min(U[idx],max(B[ichild],B[ichild+1]))
	return B
