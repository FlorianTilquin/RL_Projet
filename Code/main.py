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
from HCT_iid import HCT
from StoSOO import SOO


########## Functions to optimize ###################
def sqr(x):
    return -(x-0.5)**2

def xsin(x):
    return x*np.sin(2*np.pi*x)

def grill(x):
	if x == 0.5:
		return 1.
	u = np.abs(x-0.5)
	v = np.sqrt(u)
	s = 1.0-np.floor(2*( np.log2(u)-np.floor(np.log2(u)) ))
	return s*(v-u**2)-v+1

def graland(x):
	return 4*(1-x)*(0.75+0.25*(1-np.sqrt(np.sin(60*x))))

def sinprod(x):
	return (np.sin(13*x)*np.sin(27*x)+1.)*0.5

f = sinprod
####################################################

# Theoretical argmax
#g = lambda x:-f(x)
#Xm = fmin(g,x0 =0.8,xtol = 10**-12)
#Xm = 0.867526208796
#Ym = sinprod(Xm)
#print Xm,Ym

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
n = 5000
T,N,mu,P_HOO,REW_HOO,B = HOO2(f,1.,0.5,n)
[hm,im] = T[np.argmax(mu)]
xm = (2*im-1.)/2**(hm+1)
ym = f(xm)
print xm,ym

#L = [(2*i-1.)/2**(h+1) for [h,i] in T]
#H = [len([h for [h,i] in T if h == t]) for t in range(15)]

#print "L:",L
#print "P:",P
#print "N:",N
#print H

#### SOO ####

##n = 10000
##k = int(n/np.log(n)**3)+1
##hmax = int(np.sqrt(n/k))
##delta = np.sqrt(1./n)
##xm = SOO(g,n,k*3,hmax,delta)
##print xm

#### HCT_iid ####
# Parameters
ro = 0.5
v1 = 1.0
c1 = (ro/(3*v1))**(1.0/8)
c = 2*(ro/(1.0-ro))**(0.5)
d = 2.
Tmax = n
Ym = 1.

## Algorithm call
MC,ymc = 1,-np.inf

for i in range(MC):
	xs,ys,Xs,REWs= HCT(f,v1,ro,c,d,Tmax)
	print xs,ys
	if(ys>ymc):
		xmc=xs
		ymc=ys
		P_HCT=Xs
		REW_HCT = REWs
u = np.linspace(0,1,10001)
v = [f(x) for x in u]
plt.figure(100)
plt.plot(xmc,ymc,'bo')
plt.plot(xm,ym,'ro')
plt.plot(u,v,'-k')


plt.figure(1)
plt.plot(Ym-np.array(REW_HCT),'b+')
plt.plot(Ym-np.array(REW_HOO),'r+')

plt.figure(2)
plt.plot(P_HCT,'b+')
plt.plot(P_HOO,'r+')

V_HCT = np.arange(Tmax)*Ym-np.cumsum(REW_HCT)
V_HOO = np.arange(n+1)*Ym-np.cumsum(REW_HOO)

plt.figure(3)
plt.plot(V_HCT,'b+')
plt.plot(V_HOO,'r+')
plt.show()
