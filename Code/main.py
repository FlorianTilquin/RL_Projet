#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import numpy.random as rd
import matplotlib.pyplot as plt
from ATB import ATB
from HOO import HOO

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

## Results display
function = grill
# ATB
x = np.linspace(0,1,100)
y = function(x)

plt.figure(1)
plt.plot(x,y)
plt.plot(x[np.argmax(y)],np.max(y),'ro')
#MC,xmc,ymc = 5,0,-np.inf
#
#for i in range(MC):
#    xs,ys = ATB(function,depth,eps,gamma)
#    if(ys>ymc):
#        xmc=xs
#        ymc=ys
#plt.plot(xmc,ymc,'bo')
#plt.show()

T,N,mu = HOO(function, 1.0, 0.5, 1000)
#Ab = [(2*i-1)/2.**(h+1) for [h,i] in T]
##for x in Ab:
##	plt.plot([x,x],[-0.7,0],'r')
plt.show()
n = np.argmax(mu)
h,i = T[n]
Xm= (2*i-1.)/2**(h+1)

