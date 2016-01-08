#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import numpy.random as rd
import time
import matplotlib.pyplot as plt
from scipy.optimize import fmin
from ATB import ATB
from HOO import HOO
from HOO2 import HOO2
from POO import POO
#from StoSOO import SOO


########## Functions to optimize ###################
def sqr(x):
    return -(x-0.5)**2

def xsin(x):
    return x*np.sin(2*np.pi*x)

def grill(x):
	u = np.abs(x-0.5)
	v = np.sqrt(u)
	s = 1.0-np.floor(2*( np.log2(u)-np.floor(np.log2(u)) ))
	return s*(v-u**2)-v

f = lambda x: (np.sin(13*x)*np.sin(27*x)+1.)*0.5

####################################################

# Theoretical argmax
#g = lambda x:-f(x)
#Xm = fmin(g,x0 =0.8,xtol = 10**-12)
Xm = 0.867526208796
Ym = f(Xm)
print Xm,Ym

## Brute force
#x = np.linspace(0,1,100000)
#y = f(x)
#XM = x[np.argmax(y)]
#YM = np.max(y)
#print XM,YM

#plt.figure(1)
#plt.plot(x,y)
#plt.plot(XM,np.max(y),'ro')
#plt.show()

###################################################

######### Noise ###################################
#noise = 0.1
#g = lambda x: f(x) + rd.randn()*noise

######### Algorithms application ##################
####  ATB ####
## Parameters
#eps = 0.05
#gamma = 2.0-10**-1 ##must be smaller than 2 !!
#depth = 10
#n = 1000
#MC,xmc,ymc = 5,0,-np.inf
#
#for i in range(MC):
#    xs,ys,R = ATB(f,depth,eps,gamma,n)
#    if ys > ymc:
#        #xmc = xs
#        #ymc = ys
#		REW = R
#print xmc,ymc

####  HOO  ####
n = 1000
T,N,mu,P,REW = HOO2(f,0.,0.7,n)
[hm,im] = T[np.argmax(mu)]
xm = (2*im-1.)/2**(hm+1)
ym = f(xm)
print xm,ym
print T

##n = 10000
##k = int(n/np.log(n)**3)+1
##hmax = int(np.sqrt(n/k))
##delta = np.sqrt(1./n)
##xm = SOO(g,n,k*3,hmax,delta)
##print xm

#plt.figure(1)
#plt.plot(REW,'b+')
#plt.show()
#
#plt.figure(1)
#plt.plot(P,'b+')
#plt.show()
#
#V = np.arange(1,n+1)*Ym-np.cumsum(REW)
#
#plt.figure(3)
#plt.plot(V,'b+')
#plt.show()
