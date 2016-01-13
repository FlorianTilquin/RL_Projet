#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import numpy.random as rd

def rew(f,node) :
	H,I = node
	alpha,beta = (I-1.)/2.**H,I/2.**H
	X = rd.uniform(alpha,beta)
	return f(X),X

def SOO(f,n,k,hmax,delta,noise):
	param = np.log(n*k/delta)
	T = [[0,1]]	#Initialisation classique
	t = 0 # Nombre d'evaluataions de la fonction effectuees
	MU= [0] # La moyenne des evaluations de la fonction au noeud correspondant dans T (MU[i] = moy des Y pour T[i])
	N = [0] # Nombre d'evaluations faites au noeud T[i]
	b = [np.inf] # B values calculees sur la base de MU + coefficient correctif
	REW = []
	P = []
	depth = 0
	while t < n :
		if t%50 == 0:
			print t
		bmax = -np.inf
		for h in xrange(min(depth,hmax)+1) :
			#print "T",T
			bmax2 = -np.inf
			Lh = [j for j,[p,i] in enumerate(T) if( p == h and ([h+1,2*i] not in T) )] # tous les (indices des) feuilles de T a profondeur h
			if Lh != []:
				for i in Lh :# on calcule les b-values de tous les noeuds a profondeur h
					if b[i]>=bmax2 :
						im = i
						bmax2 = b[i]
				#print B
				if bmax2 >= bmax :
					#print "bm",b[im]
					[H,I] = T[im]
					#print im,[H,I],bmax2
					if N[im] < k :# si on a fait moins de k evaluations a ce noeud et que c'est le meilleur, on continue
						r,X = rew(f,[H,I]) # On prend un X dans la boite correspondante et on tire le bras
						r = r + rd.rand()*noise-noise/2.
						REW.append(r)
						P.append(X)
						N[im] = N[im] + 1
						MU[im] = (1-1./N[im])*MU[im] + r/N[im]
						b[im] = MU[im] + np.sqrt(param/(2.*N[im])) # calcul a la UCB
						t = t+1
						hp,ip = H-1,(I+1)/2
						while [hp,ip] in T:
							Ip = T.index([hp,ip])
							N[Ip] = N[Ip] + 1
							MU[Ip] = (1-1./N[Ip])*MU[Ip] + r/N[Ip]
							b[Ip] = MU[Ip] + np.sqrt(param/(2.*N[Ip]))
							hp,ip = hp-1,(ip+1)/2
					else :
						for s in [[H+1,2*I-1],[H+1,2*I]]:
							if s not in T:
								T.append(s)
								alpha,beta = (s[1]-1.)/2.**s[0],s[1]/2.**s[0]
								Pchild = [i for i,j in enumerate(P) if (j>=alpha)*(j<=beta)]
								if Pchild != []:
									#print "deja vu",len(Pchild),s
									N.append(len(Pchild))
									mu = np.mean([REW[i] for i in Pchild])
									MU.append(mu)
									b.append(mu  + np.sqrt(param/(2.*len(Pchild))) )
								else :
									N.append(0)
									MU.append(0)
									b.append(np.inf)
						bmax = b[im]
						if depth != H+1:
							depth = H+1
	hm = min(depth,hmax) # on va chercher la profondeur maximale atteinte pour laquelle on a calcule toutes les b-values (il se peut que depth = hmax+1 et dans ce cas les dernieres b-values n'ont pas ete calculees)
	MUh = [[MU[j],i] for j,[h,i] in enumerate(T) if h == hm] # Je vais prendre le noeud a profondeur hm qui a la plus grosse MU
	MUh = np.array(MUh)
	im = MUh[np.argmax(MUh[:,0]),1]
	#print len(T)
	return (2*im-1.)/2**(hm+1),REW,P # Milieu du meilleur segment
