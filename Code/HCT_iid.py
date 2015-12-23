import numpy as np
import numpy.random as npr
import math as m
import matplotlib.pyplot as plt
import copy

ro = 0.5
v1 = 1.0
c1 = (ro/(3*v1))**(1.0/8)
c = 2*(ro/(1.0-ro))**(0.5)
d = 1.0 ##cherche la bonne valeur
depth = 10
Tmax = 200
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
    
def sinprod(x):
    return (np.sin(13*x)*np.sin(27*x)+1.)*0.5
    
##ALGORITHM#################################################################### 
def HCT(f,depth,v1,ro,c,d,method):
    t = 1
    covTree = [[0,1],[1,1],[1,2]]
    U = np.zeros(2**depth)
    U[index(covTree[1])] = np.inf
    U[index(covTree[2])] = np.inf
    B = np.zeros(2**depth)
    T = np.zeros(2**depth)
    mu = np.zeros(2**depth)-1.0
    while(t<=Tmax):
        if(t==plus(t)):
            for I in covTree:
                idx = index(I)
                h = I[0]
                U[idx] = mu[idx]+v1*ro**h+np.sqrt(c**2*m.log(1.0/delta(t))/T[idx])
            ct = copy.copy(covTree)
            ct.reverse()
            for I in ct:
                idx = index(I)
                h,i = I[0],I[1]
                if isLeaf(covTree,[h,i]):
                    B[idx] = U[idx]
                else:
                    Bmax = max(B[index([h+1,2*i-1])],B[index([h+1,2*i])])
                    B[idx] = min(U[idx],Bmax)                
        [h,i],P = OptTraverse(covTree,B,T)
        idx = index([h,i])
        if(method=='iid'):
            left, right = interval(idx)       
            x = npr.uniform(left,right)
            r = f(x)
            t+=1
        if(method=='gamma'):
            return
        T[idx]+=1
        mu[idx] = (mu[idx]*(T[idx]-1)+r)/T[idx]
        U[idx] = mu[idx]+v1*ro**h+np.sqrt(c**2*m.log(1.0/delta(t))/T[idx])
        UpdateB(covTree,P,h,i,B,U)
        if(T[idx]>=to(h,t) and isLeaf(covTree,[h,i])):
            covTree.append([h+1,2*i-1])
            covTree.append([h+1,2*i])
            U[index([h+1,2*i-1])] = np.inf
            U[index([h+1,2*i])] = np.inf
    left,right = interval(np.argmax(mu))
    return 0.5*(left+right),np.max(mu)
    
def index(I):
    h,i = I[0],I[1]
    return 2**(h)+(i-1)

def interval(i):
    a = np.floor(m.log(i,2))
    b = 2**a
    w = i-b
    return w/(b+0.0),(w+1.0)/(b+0.0)
    
def plus(t):
    return 2**(np.floor(np.log(t))+1)
    
def OptTraverse(Tree,B,T):
    h,i,P,T[0],to = 0,1,[[0,1]],1,1
    while(T[index([h,i])]>=to and not(isLeaf(Tree,[h,i]))):
        if(B[index([h+1,2*i-1])]>=B[index([h+1,2*i])]):
            h,i = h+1,2*i-1
        else:
            h,i = h+1,2*i
        P.append([h,i])
    return [h,i],P
    
def UpdateB(Tree,P,h,i,B,U):
    idx = index([h,i])
    if isLeaf(Tree,[h,i]):
        B[idx]=U[idx]
    else:
        Bmax = max(B[index([h+1,2*i-1])],B[index([h+1,2*i])])
        B[idx] = min(U[idx],Bmax)
    Pt = copy.copy(P)
    Pt.remove([h,i])
    Pt.reverse()
    for I in Pt:
        h,i,idx = I[0],I[1],index(I) 
        Bmax = max(B[index([h+1,2*i-1])],B[index([h+1,2*i])])
        B[idx] = min(U[idx],Bmax)
    return
    
def delta(t):
    t = plus(t)
    return min(c1*d/t,1)
    
def isLeaf(Tree,I):
    h,i = I[0],I[1]
    if([h+1,2*i] in Tree or [h+1,2*i-1] in Tree):
        return False
    return True
    
def to(h,t):
    if (h==0):
        return 1.0
    return (c**2*m.log(1.0/delta(t)))/((v1*(ro**h))**2)
    
##RESULTS DISPLAY##############################################################
function = sinprod 
x = np.linspace(0,1,1000)
y = function(x)
plt.plot(x,y)
plt.plot(x[np.argmax(y)],np.max(y),'ro')
MC,xmc,ymc = 5,0.0,-np.inf

for i in range(MC):
    xs,ys = HCT(function,depth,v1,ro,c,d,'iid')
    if(ys>ymc):
        xmc=xs
        ymc=ys
plt.plot(xmc,ymc,'bo')