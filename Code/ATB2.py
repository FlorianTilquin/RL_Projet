#! /usr/bin/env/ python

import numpy as np
import numpy.random as rd


def ATB(f,Nev,depth,eps,gamma):
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
	dum = 2./gamma
	nu = 8*np.sqrt(dum*np.log2(dum))
	tau = 4./eps
	I = [np.inf]*3
	for t in xrange(1,Nev+1):
		IA = [I[D[ str(h)+str(i) ]] for [h,i] in A if h<=depth]
		print "I",IA
		print [T[D[ str(h)+str(i) ]] for [h,i] in A if h<=depth]
		i = np.argmax(IA)
		[hm,im] = A[i]
		print "Bt",[hm,im]
		alpha,beta = (im-1.)/2**hm,im/2.**hm
		x = rd.uniform(alpha,beta)
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
				N.append(1)
				Mu.append(rew)
				rho = 2.**(B[0]+1.)
				dum = 2.*np.sqrt(np.log(rho*(tau+1.)))
				R.append( dum )
				I.append(rew + (1.+2*nu)*dum)
				D[str(B[0])+str(B[1])] = len(T)-1
			else :
				print "oldB",B
				[Hb,Ib] = B
				idx = D[str(Hb)+str(Ib)]
				N[idx] += 1
				#print N[idx],Mu[idx]
				Mu[idx] = (1.-1./N[idx])*Mu[idx] + rew/N[idx]
				rho = 2.**(Hb+1.)
				dum = 2.*np.sqrt(np.log(rho*(tau+N[idx])/N[idx]))
				R[idx] = dum
				I[idx] = Mu[idx] + (1.+2*nu)*dum
			B = [B[0]-1,(1+B[1])/2]
		A_d = [a for a in A if a[0] < depth]
		for [h,i] in A_d :
			#print "Ad",A_d
			ind_p = D[str(h)+str(i)]
			for node in [[h+1,2*i-1],[h+1,2*i]]:
				if node not in T:
					T,R,I,Mu,N,D = add_child(node,T,R,I,Mu,N,D)
			ind_c0 = D[str(h+1)+str(2*i-1)]
			ind_c1 = D[str(h+1)+str(2*i)]
			r0 = R[ind_c0]
			r1 = R[ind_c1]
			m1 = 2*r0
			m2 = 2*r1
			m3 = abs(Mu[ind_c0]-Mu[ind_c1]) + r0 + r1
			if max(m1,m2,m3)>= nu*R[ind_p]:
				indA = [D[str(h)+str(i)] for [h,i] in A]
				if np.array([1. for ind in indA if N[ind] == 0]).sum() == 0:
					#print "blah"
					#print "violeur",[h,i],A
					A.remove([h,i])
					for B in [[h+1,2*i-1],[h+1,2*i]]:
						A.append(B)
						if B not in T:
							T,R,I,Mu,N,D = add_child(B,T,R,I,Mu,N,D)
			#print "A",A
			A_d = [a for a in A if a[0] < depth]
	return T,A,Mu,N,R,I,REW,POS

def add_child(node,T,R,I,Mu,N,D):
	T.append(node)
	N.append(0)
	Mu.append(0)
	R.append(np.inf)
	I.append(np.inf)
	D[str(node[0])+str(node[1])] = len(T)-1
	return T,R,I,Mu,N,D
