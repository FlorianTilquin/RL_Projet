#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import numpy.random as rd
import scipy.linalg as ln
from HOO import HOO

def POO(f,Nev,K,rhomax,numax):
	# En fait on cherche a faire aussi bien que HOO avec rho* et nu* qui sont les meilleurs parametres, sans savoir qui ils sont. On fait donc echanger les infos entre des "petits" HOO dont les parametres sont bornes par rhomax et mumax, sauf qu'en fait nu sera toujours egal a numax et qu'on s'interesse qu'a rho.
	Dmax = np.log(K)/np.log(1./rhomax)
	n = 1 # Nombre d'evaluations de la fonction 0 dans l'article mais je pense qu'ils etaient bourres vu la suite
	N = 1 # Nombre de HOO lances
	S = [[numax,rhomax]] # On va stocker les nu,rho utilises
	MU = [-np.inf] # alors la je suis pas bien sur de ce que mu doit etre : c'est la moyenne pour un HOO avec tels parametres, mais dans HOO t'as une moyenne par noeud, donc il faut faire quoi ? La moyenne sur les moyennes qu'on obtient avec HOO ? (c'est ce que je fais perso)
	while n < Nev :
		while N >= 0.5*Dmax*np.log(n/np.log(n)) and n<Nev : # j'ai un doute sur la capacite de python a se rendre compte de la valeur de n en permance, je ferai des tests c'est interessant. au passage note la pertinence du truc si n=0 ou 1...
			for i in xrange(N) : # on va lancer N HOO sur une grille exponentielle
				s = [numax,rhomax**(2.*N/(2*i+1))]# ce que j'appelle "grille exponentielle"
				S.append(s)
				T,mu = HOO(f,s[0],s[1],int(float(n)/N)) # On lance une HOO de n/N evaluations (ce qui en pratique fait pas beaucoup...)
				MU.append(np.mean(mu)) # la je sais pas trop quoi faire d'aautre que la moyenne... Le meilleur resultat sinon ?
			n = 2*n # ca c'est tres con si n= 0 au debut... Tu boucles et il se passe rien.
			N = 2*N

		for i,s in enumerate(S) :
			T,mu = HOO(f,s[0],s[1],1) #Une HOO de 1 c'est une evaluation au pif dans [0,1]. Donc c'est pas une HOO. Mais on va rien dire
			MU[i] = MU[i] + np.mean(mu)
		n = n+N
		N = N+1
		sm = S[np.argmax(MU)]
	return sm
