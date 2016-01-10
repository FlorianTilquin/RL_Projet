#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import numpy.random as rd


# Algorithm
def interval(i):
    i += 1
    a = np.floor(np.round(np.log2(i),10))
    #a = np.floor(np.log2(i))
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
         print x,d,sb
         print "box too small"

def rt(n,B,q,tau,depth):
    d = depth-np.floor(np.round(np.log2(B+1),10))
    ro = q**(d+1)
    return 2.*np.sqrt(np.log(ro*(tau+n))/n)

def ATB(f,Tmax,depth,eps,gamma,noise):
    q = 2.
    s = 2./gamma
    v = 8.*np.sqrt(s*np.log2(s))
    tau = 4./eps
    REW = []
    POS = []
    N = 2**depth-1
    mu = np.zeros(N)
    n = np.zeros(N)
    r = np.zeros(N)+0.5
    I = np.zeros(N)+np.inf
    Y = np.zeros(Tmax)
    A = [0]
    for t in xrange(Tmax):
        #print "A",A
        #B = I.tolist().index(np.max(I[A]))
        B = np.where( I == np.max(I[A]) )[0][0] #Le résultat de where est un array
        #print B
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
        #print "n[A]",n[A]
        #print "r[A]",r[A]
        while(not(np.min(n[A])==0 or len([k for k in A if W(k,mu,r,depth)>=v*r[k]])==0)):
            A1 = list(A)
            for b in (B for B in A1 if W(B,mu,r,depth)>v*r[B]):
                A.remove(b)
                A.extend([2*b+1,2*(b+1)])
    #agm = np.argmin(Y)
    xmax,ymax = POS[np.argmax(REW)],np.max(REW)
    return xmax,ymax,REW,POS

def W(B,mu,r,depth):
    W = 0.
    if(B+1 < 2**(depth-1)):
         W = max(2*r[2*B+1],2*r[2*(B+1)],abs(mu[2*B+1]-mu[2*(B+1)])+r[2*B+1]+r[2*(B+1)])
    return W

