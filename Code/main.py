#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy.random as rd
import matplotlib.pyplot as plt
#from scipy.optimize import fmin
from time import time as ti
from fonctions import *
from algos import *
#from ATB import ATB
#from HOO import HOO
#from HOO2 import HOO2
#from POO import POO
#from HCT_iid import HCT
#from StoSOO import SOO



f = sinprod
####################################################

#Theoretical argmax
#g = lambda x:-f(x)
#Xm = fmin(g,x0 =0.8,xtol = 10**-12)
if(f==sinprod):
    Xm = 0.867526208796
    Ym = sinprod(Xm)
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
noise = 0.
eps = 0.00001
gamma = 2.0-10**-1 ##must be smaller than 2 !!
depth = 15
n = 1000
xm,ym,REW_ATB,P_ATB = ATB(f,n,depth,eps,gamma,noise)
print xm,ym

####  HOO  ####
#n = 5000
#T,N,mu,P_HOO,REW_HOO,B = HOO2(f,1.,0.5,n)
#[hm,im] = T[np.argmax(mu)]
#xm = (2*im-1.)/2**(hm+1)
#ym = f(xm)
#print xm,ym

#L = [(2*i-1.)/2**(h+1) for [h,i] in T]
#H = [len([h for [h,i] in T if h == t]) for t in range(15)]

#print "L:",L
#print "P:",P
#print "N:",N
#print H

#### SOO ####

#n = 10000
#k = 20#int(n/np.log(n)**3)+1
#hmax = 20#int(np.sqrt(float(n)/k))
#delta = 0.01#np.sqrt(1./n)
#xm, REW_SOO, P_SOO = SOO(f,n,k*3,hmax,delta)
#ym = f(xm)
#print xm, ym

#### HCT_iid ####
# Parameters
#ro = 0.5
#v1 = 1.0
#c1 = (ro/(3*v1))**(1.0/8)
#c = 2*(ro/(1.0-ro))**(0.5)
#d = 2.
#Tmax = n
#Ym = 1.

## Algorithm call
#MC,ymc = 1,-np.inf
#
#for i in range(MC):
#	xs,ys,Xs,REWs= HCT(f,v1,ro,c,d,Tmax)
#	print xs,ys
#	if(ys>ymc):
#		xmc=xs
#		ymc=ys
#		P_HCT=Xs
#		REW_HCT = REWs
u = np.linspace(0,1,10001)
v = [f(x) for x in u]
plt.figure(100)
#plt.plot(xmc,ymc,'bo')
plt.plot(xm,ym,'ro')
plt.plot(u,v,'-k')



plt.figure(1)
plt.title('Positions')
#plt.plot(P_HCT,'b+')
#plt.plot(P_HOO,'r+')
#plt.plot(P_SOO,'g+')
plt.plot(P_ATB,np.arange(len(P_ATB)),'m+')

#
plt.figure(2)
plt.title('Marginal regret')
#plt.plot(Ym-np.array(REW_HCT),'b+')
#plt.plot(Ym-np.array(REW_HOO),'r+')
#plt.plot(Ym-np.array(REW_SOO),'g+')
plt.plot(Ym-np.array(REW_ATB),'m+')
#
#V_HCT = np.arange(Tmax)*Ym-np.cumsum(REW_HCT)
#V_HOO = np.arange(n+1)*Ym-np.cumsum(REW_HOO)
#V_SOO = np.arange(len(REW_SOO))*Ym-np.cumsum(REW_SOO)
V_ATB = np.arange(len(REW_ATB))*Ym-np.cumsum(REW_ATB)
#
plt.figure(3)
plt.title('Cumulative regret')
#plt.plot(V_HCT,'b+')
#plt.plot(V_HOO,'r+')
plt.plot(V_ATB,'m+')

plt.show()
