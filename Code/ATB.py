#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import numpy.random as rd


# Algorithm
def interval(i):
	i += 1
	#a = np.floor(np.round(np.log2(i),10))
	a = np.floor(np.log2(i))
	b = 2.**a
	w = i-b
	return w/b,(w+1.)/b

def parents(i,B):
	if i < B:
		print "children higher than parent"
		return np.inf
	i += 1
	B += 1
	par = [i]
	while(i/2 >= B and i > 0):
		par.append(i/2)
		i /= 2
	return np.array(par)-1

def smallest_box(x,d):
	k = np.floor(2**(d-1)*x)
	sb = (2**(d-1)+k-1).astype(int)
	if sb < 2**d:
		return sb
	else:
		print "box too small"

def rt(n,B,q,tau,depth):
	d = depth-np.floor(np.log2(B+1))
	ro = q**(d+1)
	return 2.*np.sqrt(np.log(ro*(tau+n))/n)

def ATB(f,Tmax,depth,eps,gamma,noise):
	q = 2.
	s = 2./gamma
	v = 8.*np.sqrt(s*np.log2(s))
	tau = eps/4.
	REW = []
	POS = []
	N = 2**depth-1
	mu = np.zeros(N)
	n = np.zeros(N)
	r = np.zeros(N)+np.inf
	I = np.zeros(N)+np.inf
	Y = np.zeros(Tmax)
	A = [0]
	for t in xrange(Tmax):
		#B = I.tolist().index(np.max(I[A]))
		B = np.where( I == np.max(I[A]) )[0][0] #Le résultat de where est un array
		l,h = interval(B)
		x = rd.uniform(l,h)
		y = f(x)+noise*rd.randn()
		POS.append(x)
		REW.append(y)
		sB = smallest_box(x,depth)
		par = parents(sB,B)
		# Updates
		mu[par] = (mu[par]*n[par]+y)/(n[par]+1.)
		n[par] += 1
		r[par] = rt(n[par],par,q,tau,depth)
		I = mu+(1.+2.*v)*r
		Y[t] = rt(n[B],B,q,tau,depth)
		if(np.min(n[A])>0):
			A = set_proper(A,n,mu,r,v,depth)
	#agm = np.argmin(Y)
	xmax,ymax = POS[np.argmax(REW)],np.max(REW)
	return xmax,ymax,REW,POS

def W(b,mu,r,depth):
	W = -np.inf
	if(np.floor(np.log2(b+1)) + 1 < depth):
		Lg = mu[2*b]-r[2*b]
		Ug = mu[2*b]+r[2*b]
		Lr = mu[2*b+1]-r[2*b+1]
		Ur = mu[2*b+1]+r[2*b+1]
		D = np.array([Lg-Ug,Lg-Ur,Lr-Ur,Lr-Ug])
		W = np.max(-D)##minus sign added by myself contrary to the paper
	return W

def set_proper(A,n,mu,r,v,depth):
	for b in (B for B in A if v*r[B]-W(B,mu,r,depth) <= 0):
		#print "r",r[b]
		#print W(b,mu,r,depth)
		A.remove(b)
		A.extend([2*(b+1),2*b+1])
	return A
