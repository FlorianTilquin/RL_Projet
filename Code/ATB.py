#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import numpy.random as rd


Tmax = 200
noise = 0.0
q = 1.0
gamma = 2.0-10**-1
v = 8*np.sqrt((2./gamma)*np.log2(2./gamma))
p = 1.0
depth = 10
eps = 0.05
tau = eps/4.

# Algorithm
def interval(i):
    i+=1
    a = np.floor(np.round(np.log2(i),10))
    b = 2**a
    w = i-b
    return w/(b+0.0),(w+1.0)/(b+0.0)

def parents(i,B):
    if i<B:
        print "children higher than parent"
        return np.inf
    i+=1
    B+=1
    par = [i]
    while(i/2>=B and i>0):
        par.append(i/2)
        i/=2
    return np.array(par)-1

def smallest_box(x,d):
    k = np.floor(2**(d-1)*x)
    sb = (2**(d-1)+k-1).astype(int)
    if sb<2**d:
        return sb
    else:
        print "box too small"
    
def rt(n,B):
    tau = eps/4.
    d = depth-np.floor(np.log2(B+1))
    ro = q**(p*(d+1))
    return 2*np.sqrt(np.log(ro*(tau+n))/n)

def ATB(f,depth,eps,gamma,Tmax):
    REW = []
    POS = []
    N = 2**depth-1
    mu = np.zeros(N)
    n = np.zeros(N)
    r = np.zeros(N)+np.inf
    I = np.zeros(N)+np.inf
    Y = np.zeros(Tmax)
    A = [0]
    for t in np.arange(Tmax):
        B = I.tolist().index(np.max(I[A]))
        l,h = interval(B)
        x = rd.uniform(l,h)
        y = f(x)+noise*rd.randn()
        POS.append(x)
        REW.append(y)
        sB = smallest_box(x,depth)
        par = parents(sB,B)
        # Updates
        n[par] +=1
        mu[par]=(mu[par]*(n[par]-1)+y)/n[par]
        r[par] = rt(n[par],par)
        I = mu+(1+2*p*v)*r
        Y[t] = rt(n[B],B)
        if(np.min(n[A])>0):
            A = set_proper(A,n,mu,r)
    agm = np.argmin(Y)
    xmax,ymax = POS[np.argmax(REW)],np.max(REW)
    return xmax,ymax,REW,POS

def W(b,mu,r):
    W = -np.inf
    if(np.floor(np.log2(b+1))+1<depth):
        b = b
        Lg = mu[2*b]-r[2*b]
        Ug = mu[2*b]+r[2*b]
        Lr = mu[2*b+1]-r[2*b+1]
        Ur = mu[2*b+1]+r[2*b+1]
        D = np.array([Lg-Ug,Lg-Ur,Lr-Ur,Lr-Ug])
        W=np.max(-D)##minus sign added by myself contrary to the paper ### rageux va
    return W

def set_proper(A,n,mu,r):
    for b in [B for B in A if v*r[B]-W(B,mu,r) <= 0]:
            A.remove(b)
            A.append(2*(b+1))
            A.append(2*b+1)
    return A
