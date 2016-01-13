#! /usr/bin/env python
# -*- encoding : utf-8 -*-

def plot_in_file(X,Y,filename):
	file = open(str(filename),'w')
	for x,y in zip(X,Y):
		file.write(str(x)+" "+str(y)+"\n")
	file.close()
	return
