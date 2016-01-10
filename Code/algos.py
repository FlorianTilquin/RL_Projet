#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import numpy.random as rd
from ATB import *

def HOO(f,nu,rho,Nev,T=[],N=[],mu=[],B=[],U=[],D = {},Pos =[],REW = []):
	lnev = 2*np.log(Nev)
	nuro = nu*rho**np.arange(25)
	Lnev = np.sqrt(lnev/np.arange(1,Nev+2))
	for n in xrange(1,Nev+1):
		P = [[0,1]]
		h,i = 0,1
		while [h,i] in T:
			j = D[str(h)+str(i)]
			[b1,b2] = B[j]
			[h,i] = T[j]
			if b2 > b1 :
				h,i = h+1,2*i-1
			if b2 < b1 :
				h,i = h+1,2*i
			if b2 == b1 :
				Z = rd.randint(2)
				h = h+1
				i = 2*i - Z
			P.append([h,i])
		H,I = h,i
		a,b = (I-1.)/2.**H,I/2.**H
		X = a + (b-a)*rd.rand()
		Pos.append(X)
		Y = f(X)
		REW.append(Y)
		T.append([H,I])
		D[str(H)+str(I)] = n-1
		N.append(1)
		mu.append(Y)
		U.append(mu)
		for [k,l] in P :
			n_ind = D[str(k)+str(l)]
			N[n_ind] += 1
			mu[n_ind] *= (1-1./N[n_ind])
			mu[n_ind] += Y/N[n_ind]
			U[n_ind] = mu[n_ind] + Lnev[N[n_ind]-1] + nuro[k]
		B.append( [np.inf,np.inf] )
		P1 = list(P)
		while P1 != [[0,1]] :
			[h,i] = P1[-1]
			hp,ip = h-1,(i+1)/2
			ind = D[str(hp)+str(ip)]
			B[ind][2*ip-i] = np.min([U[ind],np.max(B[ind])]) # gauche : i = 2*ip, droite: i = 2*ip-1
			P1.pop()
	return T,N,mu,B,U,D,Pos,REW

def POO(f,Nev,K,rhomax,numax):
	Dmax = np.log(K)/np.log(1./rhomax)
	n = 0 # Nombre d'evaluations de la fonction. Faut pas mettre 0 sinon ca tourne a l'infini et j'ai la flemme de faire un truc qui corrige cette connerie
	N = 1 # Nombre de HOO lances
	S = [[numax,rhomax]] # Parametres des HOO lances
	MU = [-np.inf]
	H = [[]] #Stock les HOO
	k = 0
	def condition(N,n):
		if n < 2 :
			return True
		else :
			b = 0.5*Dmax*np.log(n/np.log(n))
			return (N >= b)*(N <= 2*b)
	while n < Nev :
		print n,N
		while condition(N,n) and n < Nev :
			for i in xrange(N) :
				s = [numax,rhomax**(2.*N/(2*i+1))]
				s = [numax,0.5]
				S.append(s)
				NevH = int(float(n)/N)+1
				dum = list(HOO(f,s[0],s[1],NevH))
				H.append( dum ) # On lance une HOO de n/N evaluations (ce qui en pratique fait pas beaucoup...)
				MU.append(np.max(dum[2])) # On prend le meilleur resultat de cette HOO
			n = n + N * NevH
			N = 2*N
			k = k + NevH
		for i,s in enumerate(S) :
			H[i] = list(HOO(f,s[0],s[1],1,*H[i]))
			#Ni = H[i][1]
			MU[i] = np.max(H[i][2]) #MU[i]*(1-1./k) + np.max(H[i][2])/k
		n = n+N
		k = k+1
		sm = np.argmax(MU)
		#print "s,S[s]:",sm,S[sm]
	return H#list(H[sm])

#################################################### SOO ###########################################

def rew(f,node) :
    H,I = node
    alpha,beta = (I-1.)/2.**H,I/2.**H
    X = alpha + (beta-alpha)*rd.rand()
    rew = f(X) # On prend un X dans la boite correspondante et on tire le bras
    return rew,X

def SOO(f,n,k,hmax,delta):
    T = [[0,1]]    #Initialisation classique
    t = 0 # Nombre d'evaluataions de la fonction effectuees
    MU= [0] # La moyenne des evaluations de la fonction au noeud correspondant dans T (MU[i] = moy des Y pour T[i])
    N = [0] # Nombre d'evaluations faites au noeud T[i]
    b = [np.inf] # B values calculees sur la base de MU + coefficient correctif
    REW = []
    P = []
    while t < n :
        #print "enter main loop"
        #print "t=",t
        bmax = -np.inf
        depth = np.max([h for [h,i] in T])
        for h in xrange(min(depth,hmax)+1) :
            #print "T",T
            #print "enter B-loop"
            B = []
            Lh = [j for j,[p,i] in enumerate(T) if p == h and 1-([h+1,2*i] in T)] # tous les (indices des) noeuds de T a profondeur h
            if Lh != []:
                #print T
                for i in Lh :# on calcule les b-values de tous les noeuds a profondeur h
                    if N[i] == 0 :
                        #print "new node initialized at t=",t
                        b[i] = np.inf # +inf si on y est jamais passe, pour forcer au moins un passage
                    else :
                        b[i] = MU[i] + np.sqrt(np.log(n*k/delta)/(2*N[i])) # calcul a la UCB
                    B.append(b[i])
                #print B
                im = Lh[np.argmax(B)] # Index of the best node in T at depth h
                #print h,im
                if b[im] >= bmax :
                    [H,I] = T[im]
                    if N[im] < k :# si on a fait moins de k evaluations a ce noeud et que c'est le meilleur, on continue
                        r,X = rew(f,[H,I]) # On prend un X dans la boite correspondante et on tire le bras
                        REW.append(r)
                        P.append(X)
                        N[im] = N[im] + 1
                        MU[im] = (1-1./N[im])*MU[im] + r/N[im]
                        t = t+1
                    else :
                        for s in [[H+1,2*I],[H+1,2*I-1]]:
                            if s not in T:
                                T.append(s)
                                b.append(np.inf)
                                MU.append(0)
                                N.append(0)
                            else :
                                #MU[im] = -np.inf
                                r,X = rew(f,[H,I]) # On prend un X dans la boite correspondante et on tire le bras
                                REW.append(r)
                                P.append(X)
                                N[im] = N[im] + 1
                                MU[im] = (1-1./N[im])*MU[im] + r/N[im]
                                t = t+1
                    bmax = b[im]
    hm = min(depth,hmax) # on va chercher la profondeur maximale atteinte pour laquelle on a calcule toutes les b-values (il se peut que depth = hmax+1 et dans ce cas les dernieres b-values n'ont pas ete calculees)
    MUh = [[MU[j],i] for j,[h,i] in enumerate(T) if h == hm] # Je vais prendre le noeud a profondeur hm qui a la plus grosse MU
    MUh = np.array(MUh)
    im = MUh[np.argmax(MUh[:,0]),1] # et voila je l'ai fait
    return (2*im-1.)/2**(hm+1),REW,P # Milieu du meilleur segment


############# HCT #########################################################
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
				U[idx][1] = U[idx][-1] + v1*ro**h + c*np.sqrt(np.log(1.0/delta[t-1])/U[idx][0])
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
