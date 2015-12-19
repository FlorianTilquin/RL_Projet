#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import numpy.random as rd
import matplotlib.pyplot as plt

##PARAMETERS###################################################################
eps = 0.05
gamma = 2.0-10**-1 ##must be smaller than 2 !!
p = 1.0
v = 8*np.sqrt((2./gamma)*np.log2(2./gamma))
tau = 4./eps
q = 1.0
depth = 10
Tmax = 200
noise = 0.0

##FUNCTIONS TO OPTIMIZE########################################################
def sqr(x):
    return -(x-0.5)**2

def xsin(x):
    return x*np.sin(2*np.pi*x)

def grill(x):
	u = np.abs(x-0.5)
	v = np.sqrt(u)
	s = 1.0-np.floor(2*( np.log2(u)-np.floor(np.log2(u)) ))
	return s*(v-u**2)-v

##ALGORITHM####################################################################
def dyadiqueTree(depth):
    return np.arange(2**depth-1)+1

def interval(i):
    a = np.floor(np.log2(i))
    b = 2**a
    w = i-b
    return w/(b+0.0),(w+1.0)/(b+0.0)

def parents(i,B):
    K = np.floor(np.log2(i))+1
    L = np.floor(np.log2(B))+1
    prof = K-L+1
    Den = (np.zeros(prof)+2)**(np.arange(prof))
    return np.floor((np.zeros(prof)+i)/Den).astype(int)-1

def smallest_box(x,Tree):
    N = np.floor(np.log2(Tree[-1]))+1
    k = np.ceil(2**(N-1)*x)
    return Tree[2**(N-1)+k-2]

def ATB(f,depth,eps,gamma):
	Tree = dyadiqueTree(depth)
	#print "Tree : ",Tree
	N = len(Tree)
	mu = np.zeros(N)
	n = np.zeros(N)
	r = np.zeros(N)+np.inf
	I = np.zeros(N)+np.inf
	Y = np.zeros([Tmax,3])
	A = [1]
	for t in np.arange(Tmax):
		B = A[np.argmax(I[np.array(A)-1])]
		l,h = interval(B)
		x = rd.uniform(l,h)
		y = f(x)+noise*rd.randn()
		sB = smallest_box(x,Tree)
		par = parents(sB,B)
		# Updates
		n[par] +=1
		mu[par]=(mu[par]*(n[par]-1)+y)/n[par]
		r[par] = rt(n[par],par)
		I = mu+(1+2*p*v)*r
		Y[t,0],Y[t,1] = x,y
		Y[t,2] = rt(n[B-1],B)
		if(np.min(n[np.array(A)-1])>0):
			A = set_proper(A,n,mu,r)
	agm = np.argmin(Y[:,2])
	xmax,ymax = Y[agm,0],Y[agm,1]
	return xmax,ymax

def rt(n,B):
	d = depth-np.floor(np.log2(B))
	ro = q**(p*(d+1))
	return 2*np.sqrt(np.log(ro*(tau+n))/n)

def W(B,mu,r):
    if(type(B)==int):
        B = [B]
        W = np.zeros(1)
    else:
        W = np.zeros(len(B))
    t=0
    for b in B:
        if(np.floor(np.log2(b))+1<depth):
            b = b-1
            Lg = mu[2*b]-r[2*b]
            Ug = mu[2*b]+r[2*b]
            Lr = mu[2*b+1]-r[2*b+1]
            Ur = mu[2*b+1]+r[2*b+1]
            D = np.array([Lg-Ug,Lg-Ur,Lr-Ur,Lr-Ug])
            W[t]=np.max(-D)##minus sign added by myself contrary to the paper ### rageux va
        t+=1
    return W

def set_proper(A,n,mu,r):
    for b in A:
		if np.min(v*r[b-1]-W(b,mu,r)) <= 0 :
			A.remove(b)
			A.append(2*b)
			A.append(2*b+1)
    return A

##RESULTS DISPLAY##############################################################
function = grill
x = np.linspace(0,1,100)
y = function(x)

plt.figure(1)
plt.plot(x,y)
plt.plot(x[np.argmax(y)],np.max(y),'ro')
MC,xmc,ymc = 5,0,-np.inf

for i in range(MC):
    xs,ys = ATB(function,depth,eps,gamma)
    if(ys>ymc):
        xmc=xs
        ymc=ys
plt.plot(xmc,ymc,'bo')
plt.show()
