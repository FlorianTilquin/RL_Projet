#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy.random as rd
import matplotlib.pyplot as plt
from scipy.optimize import fmin
from time import time as ti
from fonctions import *
#from algos import *
from ATB2 import ATB
from HOO import HOO
from HCT2 import *
#from HOO2 import HOO2
from POO import POO
#from HCT_iid import HCT
#from StoSOO import SOO



#f = sinprod
#f = sqr
f = sinprod

#Theoretical argmax
#g = lambda x:-f(x)
#Xm = fmin(g,x0 =0.5235,xtol = 10**-12)
#print Xm
if(f==sinprod):
	Xm = 0.867526208796
	Ym = f(Xm)
	print Xm,Ym

#if(f==garland):
#	#Xm = 0.523598784402
#	Xm = 0.523598775599
#	Ym = f(Xm)
#	print Xm,Ym

# Brute force
#x = np.linspace(0,1,10000)
#y = f(x)
#Xm = x[np.argmax(y)]
#Ym = np.max(y)
#print Xm,Ym

#plt.figure(1)
#plt.plot(x,y)
#plt.plot(Xm,np.max(y),'ro')
#plt.show()

###################################################

######### Noise ###################################
#noise = 0.01
#f = lambda x: g(x) + rd.randn()*noise


######### Algorithms application ##################
####  ATB ####
## Parameters
noise = 0.
eps = 0.01
gamma = 1.99
depth = 15
n = 10
T,A,Mu,N,R,I,REW_ATB,P_ATB = ATB(f,n,depth,eps,gamma)#,noise)
#print len(T),len(A)
[hm,im] = T[np.argmax(REW_ATB)]
xm = (2*im-1.)/2**(hm+1)
ym = f(xm)
print xm,ym

####  HOO  ####
#n = 5000
#T,N,mu,B,U,D,P_HOO,REW_HOO = HOO(f,1.,0.5,n)
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

#### POO ####
#n = 5000
#K = 2
#rhomax = 0.8
#numax = 1.
#H = POO(f,n,K,rhomax,numax)
#T = H[0]
#mu = H[2]
#print len(mu)
#P_POO = H[6]
#REW_POO = H[7]
#im = np.argmax(mu)
#[h,i] = T[im]
#xm = (2*i-1.)/2**(h+1)
#ym = f(xm)
#print Xm,Ym
#print xm,ym

#### SOO ####

#n = 5000
#k = (int(n/np.log(n)**3))
#hmax = int(np.log(n)**1.5)
#delta = np.sqrt(1./n)
#xm, REW_SOO, P_SOO = SOO(f,n,k,hmax,delta)
#ym = f(xm)
#print xm, ym
#print Xm,Ym

#### HCT_iid ####
# Parameters
#n = 1000
#rho = 0.7
#nu = 1.0
#delta = 0.1
#
#T,N,Mu,U,B,D,P_HCT,REW_HCT = HCT(f,n,nu,rho,delta)
#print N
#print len(T)
#im = np.argmax(Mu)
#[h,i] = T[im]
#xm = (2*i-1.)/2**(h+1)
#ym = f(xm)
#print xm,ym

u = np.linspace(0,1,10001)
v = [f(x) for x in u]
plt.figure(100)
plt.plot(xm,ym,'ro')
plt.plot(Xm,Ym,'bo')
plt.plot(u,v,'-k')



plt.figure(1)
plt.title('Positions')
#plt.plot(P_HCT,np.arange(len(P_HCT)),'b+')
#plt.plot(P_HOO,np.arange(len(P_HOO)),'r+')
#plt.plot(P_POO,np.arange(len(P_POO)),'c+')
#plt.plot(P_SOO,np.arange(len(P_SOO)),'g+')
plt.plot(P_ATB,np.arange(len(P_ATB)),'m+')

#Ym = 0.
plt.figure(2)
plt.title('Marginal regret')
#plt.plot(np.array(REW_HCT),'b+')
#plt.plot(Ym-np.array(REW_HOO),'r+')
#plt.plot(Ym-np.array(REW_POO),'c+')
#plt.plot(Ym-np.array(REW_SOO),'g+')
plt.plot(Ym-np.array(REW_ATB),'m+')
#
#V_HCT = np.arange(n)*Ym-np.cumsum(REW_HCT)
##V_HOO = np.arange(n)*Ym-np.cumsum(REW_HOO)
##V_SOO = np.arange(n)*Ym-np.cumsum(REW_SOO)
#V_POO = np.arange(len(REW_POO))*Ym-np.cumsum(REW_POO)
print Ym
V_ATB = np.arange(n)*Ym-np.cumsum(REW_ATB)
##
plt.figure(3)
plt.title('Cumulative regret')
#plt.plot(V_HCT,'b+')
##plt.plot(V_HOO,'r+')
#plt.plot(V_POO,'c+')
##plt.plot(V_SOO,'r+')
plt.plot(V_ATB,'m+')
#
plt.show()
