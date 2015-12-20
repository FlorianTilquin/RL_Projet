#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import numpy.random as rd
from time import time as ti
import matplotlib.pyplot as plt
#from ATB import ATB
from HOO import HOO
from POO import POO

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
#function = lambda x: (np.sin(13*x)*np.sin(27*x)+1.)*0.5
function = grill
# ATB
# x = np.linspace(0,1,10000)
# y = function(x)
# 
# plt.figure(1)
# plt.plot(x,y)
# plt.plot(x[np.argmax(y)],np.max(y),'ro')
#MC,xmc,ymc = 5,0,-np.inf
#
#for i in range(MC):
#    xs,ys = ATB(function,depth,eps,gamma)
#    if(ys>ymc):
#        xmc=xs
#        ymc=ys
#plt.plot(xmc,ymc,'bo')
#plt.show()

sm = POO(function,1000, 2,	0.9, 1.2)
print sm
T,mu = HOO(function, 1.0, 0.5, 100)

#Ab = [(2*i-1)/2.**(h+1) for [h,i] in T]
#for x in Ab:
#	plt.plot([x,x],[-0.7,0],'r')
# plt.show()
# n = np.argmax(mu)
# h,i = T[n]
# Xm = (2*i-1.)/2**(h+1)
# plt.plot(Xm,function(Xm),'bo')
