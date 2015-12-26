#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import numpy.random as rd
import time
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
print(x[np.argmax(y)])
noise = 0.1
g = lambda x: f(x) + rd.randn()*noise

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

# n = np.argmax(mu)
# h,i = T[n]
# Xm = (2*i-1.)/2**(h+1)
# plt.plot(Xm,f(Xm),'bo')

n = 10000
k = int(n/np.log(n)**3)+1
hmax = int(np.sqrt(n/k))
delta = np.sqrt(1./n)
xm = SOO(g,n,k*3,hmax,delta)
print xm
## POO
#a = time.time()
#T,N,mu,B,U,D = POO(f,n,2,0.9,1.)
#b = time.time()
#print b-a
#[hm,im] = T[np.argmax(mu)]
#print (2*im-1.)/2**(hm+1)
