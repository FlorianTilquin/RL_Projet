# -*- coding: utf-8 -*-
"""
Created on Wed Dec 09 15:17:30 2015

@author: Nathan
"""
import numpy as np
import numpy.random as npr
import math as m
import matplotlib.pyplot as plt

##PARAMETERS###################################################################
eps = 0.05
gamma = 2.0-10**(-1) ##must be smaller than 2 !!
p = 1.0
v = 8*np.sqrt((2/gamma)*m.log((2/gamma),2))
tau = 4/eps
q = 1.0
depth = 10
Tmax = 200
noise = 0.0
##FUNCTIONS TO OPTIMIZE########################################################
def sqr(x):
    return -(x-0.5)**2
    
def xsin(x):
    return x*np.sin(2*np.pi*x)
    
def grill(x):
    u = x-0.5
    def s(v):
        w = v-np.floor(v)
        return 1.0-np.floor(2*w)
    return s(np.log(abs(u))/np.log(2))*(np.sqrt(abs(u))-u**2)-np.sqrt(abs(u))
    
##ALGORITHM####################################################################    
def dyadiqueTree(depth):
    return np.arange(2**depth-1)+1
    
def interval(i):
    a = np.floor(m.log(i,2))
    b = 2**a
    w = i-b
    return w/(b+0.0),(w+1.0)/(b+0.0)
    
def parents(i,B):
    K = np.floor(m.log(i,2))+1
    L = np.floor(m.log(B,2))+1
    prof = K-L+1
    Den = (np.zeros(prof)+2)**(np.arange(prof))
    return np.floor((np.zeros(prof)+i)/Den).astype(int)
    
def smallest_box(x,Tree):
    N = np.floor(m.log(Tree[-1],2))+1
    k = np.ceil(2**(N-1)*x)
    return Tree[2**(N-1)+k-2]

def ATB(f,depth,eps,gamma):
    Tree = dyadiqueTree(depth)
    #print "Tree : ",Tree
    N = len(Tree)
    mu = np.zeros(N)
    n = np.zeros(N)
    r = np.zeros(N)+np.inf
    I = np.zeros(N)+np.inf
    Y = np.zeros([Tmax,3])        
    A = [1]
    for t in np.arange(Tmax):
        #print 
        #print "active boxes :",A
        #print "It|A: ",I[np.array(A)-1]
        B = A[np.argmax(I[np.array(A)-1])]
        #print "B = ",B
        l,h = interval(B)
        x = npr.uniform(l,h)
        #print "arm pulled: ",x
        y = f(x)+noise*npr.randn()        
        sB = smallest_box(x,Tree)
        #print "smallest box corresponding: ",sB
        par = parents(sB,B)
        #print "parents of the smallest box: ",par
        n = update_n(n,par)
        #print "n =",n
        if(n[B-1]==0):
            return "Error in the counts"
        mu = update_mu(mu,n,par,y)
        #print "mu =",mu
        r = update_r(r,n,par)
        #print "r =",r
        I = It(mu,r)
        if(I[B-1]==np.inf):
            return "Error"
        #print "I =",I
        Y[t,0],Y[t,1],Y[t,2] = x,y,rt(n[B-1],B)
        if(np.min(n[np.array(A)-1])>0):
            #print "check if A is proper"
            A = set_proper(A,n,mu,r)
    agm = np.argmin(Y[:,2])
    xmax,ymax = Y[agm,0],Y[agm,1]               
    return xmax,ymax
    
def update_n(n,P):
    n[P-1]+=1
    return n
    
def update_mu(mu,n,P,y):
    P = P-1
    mu[P]=(mu[P]*(n[P]-1)+y)/n[P]
    return mu
    
def update_r(r,n,P):
    r[P-1]=rt(n[P-1],P)
    return r
   
def It(mu,r):
    return mu+(1+2*p*v)*r
    
def rt(n,B):
    return 2*np.sqrt(np.log(ro(B)*(tau+n))/n)
    
def ro(B):
    return q**(p*(d(B)+1))
    
def d(B):
    return depth-np.floor(np.log(B)/np.log(2))
    
def W(B,mu,r):
    if(type(B)==int):
        B = [B]
        W = np.zeros(1)
    else:
        W = np.zeros(len(B))
    t=0
    for b in B:        
        #print "box to be verified ",b
        #print "depth of this box :",np.floor(m.log(b,2))+1
        if(np.floor(m.log(b,2))+1<depth):
            b = b-1
            Lg = mu[2*b]-r[2*b]
            Ug = mu[2*b]+r[2*b]
            Lr = mu[2*b+1]-r[2*b+1]
            Ur = mu[2*b+1]+r[2*b+1]
            D = np.array([Lg-Ug,Lg-Ur,Lr-Ur,Lr-Ug])
            W[t]=np.max(-D)##minus sign added by myself contrary to the paper
        t+=1
    return W
    
def invariant(A,mu,r):
    ii = np.min(v*r[A-1]-W(A,mu,r))
    #print "second condition :",ii
    if(ii>=0):
        return True
    else:    
        return False
        
def set_proper(A,n,mu,r):
    for b in A:
        if(invariant(b,mu,r)==False):
            #print "modifying A"
            A.remove(b)
            A.append(2*b)
            A.append(2*b+1)
    return A

##RESULTS DISPLAY##############################################################
function = grill 
x = np.linspace(0,1,100)
y = function(x)
plt.plot(x,y)
plt.plot(x[np.argmax(y)],np.max(y),'ro')
MC,xmc,ymc = 5,0,-np.inf

for i in range(MC):
    xs,ys = ATB(function,depth,eps,gamma)
    if(ys>ymc):
        xmc=xs
        ymc=ys
plt.plot(xmc,ymc,'bo')
    

    
    

    

    

    
    
