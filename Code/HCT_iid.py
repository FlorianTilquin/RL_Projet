#! /usr/bin/env python

import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt

##ALGORITHM####################################################################
def HCT(f,v1,ro,c,d,Tmax):
	REW = []
	X = []
	par = [c,v1,ro]
	c1 = (ro/(3*v1))**(1.0/8)
	T = np.arange(Tmax+1)+1
	T = 2**(np.floor(np.log(T))+1)
	T = c1*d/T
	T = np.vstack(( T, np.ones(Tmax+1) ))
	delta = np.min( T, 0)
	t = 1
	covTree = [[0,1],[1,1],[1,2]]
	U = [[0,np.inf,np.inf,-1.0],[0,np.inf,np.inf,-1.0],[0,np.inf,np.inf,-1.0]]
	while(t<=Tmax):
		if( t==2**(np.floor(np.log(t))+1) ):
			for i in xrange(len(covTree)):
				h = covTree[i][0]
				U[i][1] = U[i][-1] + v1*ro**h + c*np.sqrt(np.log(1.0/delta[t-1])/U[i][0])
			ct = list(covTree)
			ct.reverse()
			for k in xrange(len(covTree)):
				idx = -(k+1)
				V = U[idx]
				I = covTree[idx]
				h,i = I[0],I[1]
				if isLeaf(covTree,[h,i]):
					V[2]=V[1]
				else:
					Bmax = max(U[covTree.index([h+1,2*i-1])][2],U[covTree.index([h+1,2*i])][2])
					V[2]=min(V[1],Bmax)
		[h,i],P = OptTraverse(covTree,U,t,par,delta)
		idx = covTree.index([h,i])
		left, right = (i-1)/(2.**h),i/(2.**h)
		#print left,right
		x = npr.uniform(left,right)
		X.append(x)
		r = f(x)
		REW.append(r)
		#print r
		t+=1
		V = U[idx]
		V[0]+=1
		V[-1] = (V[-1]*(V[0]-1)+r)/V[0]
		V[1] = V[-1]+v1*ro**h+c*np.sqrt(np.log(1.0/delta[t-1])/V[0])#max(V[0],1))
		UpdateB(covTree,P,h,i,U)
		to = (h==0)*1.+(h>0)*(c**2*np.log(1.0/delta[t-1]))/((v1*(ro**h))**2)
		if(V[0]>=to and isLeaf(covTree,[h,i])):# and index([h+1,2*i])<2**depth):
			covTree.append([h+1,2*i-1])
			covTree.append([h+1,2*i])
			U.append([0,np.inf,np.inf,-1.0])
			U.append([0,np.inf,np.inf,-1.0])
	U = np.array(U).T
	idx = np.argmax(U[-1])
	[h,i] = covTree[idx]
	left,right = (i-1)/(2.**h),i/(2.**h)
	#print left,right
	return 0.5*(left+right),f(0.5*(left+right)),X,REW

def OptTraverse(Tree,U,t,par,delta):
	[c,v1,ro] = par
	h,i,P = 0,1,[[0,1]]
	to = (h==0)*1.+(h>0)*(c**2*np.log(1.0/delta[t-1]))/((v1*(ro**h))**2)

	while(U[Tree.index([h,i])][0]>= to and not(isLeaf(Tree,[h,i]))):
		if(U[Tree.index([h+1,2*i-1])][2]>U[Tree.index([h+1,2*i])][2]):
			h,i = h+1,2*i-1
		elif(U[Tree.index([h+1,2*i-1])][2]<U[Tree.index([h+1,2*i])][2]):
			h,i = h+1,2*i
		else:
			Z = npr.randint(2)
			h,i = h+1,2*i-Z
		P.append([h,i])
	return [h,i],P

def UpdateB(Tree,P,h,i,U):
    idx = Tree.index([h,i])
    if isLeaf(Tree,[h,i]):
        U[idx][2]=U[idx][1]
    else:
        Bmax = max(U[Tree.index([h+1,2*i-1])][2],U[Tree.index([h+1,2*i])][2])
        U[idx][2]=min(U[idx][1],Bmax)
	Pt = list(P)
	Pt.remove([h,i])
	for k in xrange(len(Pt)):
		I = Pt[-(k+1)]
		h,i = I[0],I[1]
		idx = Tree.index(I)
        Bmax = max(U[Tree.index([h+1,2*i-1])][2],U[Tree.index([h+1,2*i])][2])
        U[idx][2]=min(U[idx][1],Bmax)
    return


def isLeaf(Tree,I):
	h,i = I[0],I[1]
	return 1 - max([h+1,2*i] in Tree , [h+1,2*i-1] in Tree)
