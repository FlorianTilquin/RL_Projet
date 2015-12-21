#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import numpy.random as rd
from time import time as ti
import matplotlib.pyplot as plt
#from ATB import ATB
from HOO import HOO
from POO import POO
from StoSOO import SOO

# Parameters
eps = 0.05
gamma = 2.0-10**-1 ##must be smaller than 2 !!
depth = 10

# Functions to optimize
def sqr(x):
    return -(x-0.5)**2

def xsin(x):
    return x*np.sin(2*np.pi*x)

def grill(x):
	u = np.abs(x-0.5)
	v = np.sqrt(u)
	s = 1.0-np.floor(2*( np.log2(u)-np.floor(np.log2(u)) ))
	return s*(v-u**2)-v

# Results display
f = lambda x: (np.sin(13*x)*np.sin(27*x)+1.)*0.5
#f = grill
# ATB
x = np.linspace(0,1,10000)
y = f(x)
print x[np.argmax(y)]

#plt.figure(1)
#plt.plot(x,y)
#plt.plot(x[np.argmax(y)],np.max(y),'ro')
#plt.show()
#MC,xmc,ymc = 5,0,-np.inf
#
#for i in range(MC):
#    xs,ys = ATB(f,depth,eps,gamma)
#    if(ys>ymc):
#        xmc=xs
#        ymc=ys
#plt.plot(xmc,ymc,'bo')
#plt.show()

#sm = POO(f,1000, 2,	0.9, 1.2)
#T,mu = HOO(f, 1.0, 0.5, 100)
# n = np.argmax(mu)
# h,i = T[n]
# Xm = (2*i-1.)/2**(h+1)
# plt.plot(Xm,f(Xm),'bo')

n = 1000.
k = int(n/np.log(n)**3)+1
hmax = int(np.sqrt(n/k))
delta = np.sqrt(1./n)
print k, hmax, delta
xm = SOO(f,n,k+50,hmax,delta)
print xm

