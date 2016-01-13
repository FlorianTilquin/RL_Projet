#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy.random as rd
import matplotlib.pyplot as plt
from scipy.optimize import fmin
from time import time as ti
from fonctions import *
from ATB2 import ATB
from HOO import HOO
from HCT2 import *
from POO import *
from StoSOO import *
from plot_in_file import *



f = sinprod
#f = sqr
#f = garland
#f = grill

#Theoretical argmax
#g = lambda x:-f(x)
#Xm = fmin(g,x0 =0.5235,xtol = 10**-12)
#print Xm
if(f==sinprod):
	Xm = 0.867526208796
	Ym = f(Xm)

if(f==garland):
	Xm = 0.523598775599
	Ym = f(Xm)

if(f== grill):
	Xm = 0.5
	Ym = 1.

######### Noise ###################################
noise = 0.
n = 1000
#f = lambda x: g(x) + rd.randn()*noise


######### Algorithms application ##################
####  ATB ####
# Default parameters
eps = 10**-6
nu = 3.2 # 2.5<nu<3.5
depth = 15
def_ATB = [depth,eps,nu,noise]

####  HOO  ####
# Default parameters
nu = 1.
rho= 0.5
def_HOO = [nu,rho,noise]

#### POO ####
rhomax = 0.8
numax = 1.
def_POO = [rhomax,numax,noise]

#### SOO ####
k = (np.ceil(n/np.log(n)**3))
hmax = np.ceil(np.log(n)**1.5)
delta = np.sqrt(1./n)
def_SOO = [k,hmax,delta,noise]

#### HCT_iid ####
# Default parameters
rho = 0.1
nu = 1.0
delta = 0.8
def_HCT = [nu,rho,delta,noise]


#for algo in [HOO]:#,POO,HCT,SOO]:
#print algo.__name__
algo = HOO
Result = algo(f,n,*eval('def_'+algo.__name__))
P,R = Result[-2],Result[-1]

plt.figure(1)
plt.title('Positions')
plt.plot(P,np.arange(len(P)),'k.')

plt.figure(2)
plt.title('Marginal regret')
plt.plot(Ym-np.array(R),'k.')

V = np.arange(len(R))*Ym-np.cumsum(R)

plt.figure(3)
plt.title('Cumulative regret')
plt.plot(V,'k.')
plt.show()

plt.figure(4)
plt.title('Mean regret')
plt.plot(V/np.arange(1,len(V)+1),'k.')
plt.show()
