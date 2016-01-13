#! /usr/bin/env/ python

import numpy as np
import numpy.random as rd


def ATB(f,Nev,depth,eps,nu):
	POS = []
	REW = []
	A = [[0,1]]
	T = [[0,1],[1,1],[1,2]]
	D = {}
	D['01'] = 0
	D['11'] = 1
	D['12'] = 2
	N = [0]*3
	Mu = [0.]*3
	R = [np.inf]*3
	#dum_rad = np.sqrt(8*np.log(Nev))
	#dum = 2./gamma
	#nu = 8*np.sqrt(dum*np.log2(dum))
	tau = 4./eps
	I = [np.inf]*3
	for t in xrange(1,Nev+1):
		print len(A)
		print "etape",t
		IA = [I[D[ str(h)+str(i) ]] for [h,i] in A if h<=depth]
		#print "I",IA
		#print "A-IA",[T[D[ str(h)+str(i) ]] for [h,i] in A if h<=depth]
		#print "n-IA",[N[D[ str(h)+str(i) ]] for [h,i] in A if h<=depth]
		i = np.argmax(IA)
		[hm,im] = A[i]
		#print "Bt",[hm,im]
		alpha,beta = (im-1.)/2**hm,im/2.**hm
		x = rd.uniform(alpha,beta)
		#print x
		POS.append(x)
		rew = f(x)
		REW.append(rew)
		ix = np.ceil(x*2**depth).astype(int)
		B = [depth,ix]
		#print "smox",B
		while B[0] >= hm:
			if B not in T:
				#print "newB",B
				T.append(B)
				#print "B j'ajoute le node",B,"a T"
				N.append(1)
				Mu.append(rew)
				rho = 2.**-(B[0]+1.)## + ou - ??
				dum = 2.*np.sqrt(np.log(rho*(tau+1.)))
				R.append( dum )
				I.append(rew + (1.+2*nu)*dum)
				D[str(B[0])+str(B[1])] = len(T)-1
			else :
				#print "oldB",B
				[Hb,Ib] = B
				idx = D[str(Hb)+str(Ib)]
				N[idx] += 1
				Mu[idx] = (1.-1./N[idx])*Mu[idx] + rew/N[idx]
				rho = 2.**-(Hb+1.)
				dum = 2.*np.sqrt(np.log( rho*(tau+N[idx]) )/N[idx])
				R[idx] = dum
				#print "N,R,Mu",N[idx],(1.+2*nu)*dum*0.01,Mu[idx]
				#print T[idx]
				I[idx] = Mu[idx] + (1.+2*nu)*dum
			B = [B[0]-1,(1+B[1])/2]
			condition = 1
		while condition == 1:
			#print "while"
			indA = [D[str(h)+str(i)] for [h,i] in A] # h< depth ?
			if np.array([1. for ind in indA if N[ind] == 0]).sum() == 0: # si tout le monde a ete visite au moins une fois
				#print "la condition 1 est valide"
				A_d = [a for a in A if a[0] < depth]
				condition = 0
				for [h,i] in A_d : # On construit le W de chacun dans A
					#print "Ad",A_d
					#print "h,i in A",[h,i],A_d
					ind_p = D[str(h)+str(i)]
					for node in [[h+1,2*i-1],[h+1,2*i]]:
						if node not in T:
							T,R,I,Mu,N,D = add_child(node,T,R,I,Mu,N,D)
							#print "A j'ajoute l'enfant ",node, "a T"
						#else :
							#print node,"est deja dans T"
					ind_c0 = D[str(h+1)+str(2*i-1)]
					ind_c1 = D[str(h+1)+str(2*i)]
					r0 = R[ind_c0]
					r1 = R[ind_c1]
					m1 = 2*r0
					m2 = 2*r1
					m3 = abs(Mu[ind_c0]-Mu[ind_c1]) + r0 + r1
					#print m3,abs(Mu[ind_c0]-Mu[ind_c1]),r0,r1
					if max(m1,m2) == np.inf:
						cond2 = min(m1,m2)
					else :
						cond2 = max(m1,m2,m3)
					#print cond2,nu*R[ind_p] # ce max est-il souvent fini dans les faits ?
					if cond2 >= nu*R[ind_p]:
						condition = 1
						#print "la conditon c2 est validee par:"
						#print "le violeur",[h,i],"dans",A
						A.remove([h,i])
						for B in [[h+1,2*i-1],[h+1,2*i]]:
							#print "enfant issus du viol",B,"ajoute a A"
							A.append(B)
							if B not in T:
								#print "Ad j'ajoute l'enfant",B,"a T"
								T,R,I,Mu,N,D = add_child(B,T,R,I,Mu,N,D)
					#else :
					#	condition = 0
					#	break
			else :
				condition = 0
				#print "Il y a un noeud qui necessite une visite; Sortie de boucle"
				break
		#	print "A",A
	return T,A,Mu,N,R,I,REW,POS

def add_child(node,T,R,I,Mu,N,D):
	T.append(node)
	N.append(0)
	Mu.append(0)
	R.append(np.inf)
	I.append(np.inf)
	D[str(node[0])+str(node[1])] = len(T)-1
	return T,R,I,Mu,N,D
