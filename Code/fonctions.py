import numpy as np

########## Functions to optimize ###################
def sqr(x):
    return -(x-0.5)**2

def xsin(x):
    return x*np.sin(2*np.pi*x)

def grill(x):
	if x == 0.5:
		return 1.
	u = np.abs(x-0.5)
	v = np.sqrt(u)
	s = 1.0-np.floor(2*( np.log2(u)-np.floor(np.log2(u)) ))
	return s*(v-u**2)-v+1

def garland(x):
	return 4*x*(1-x)*(0.75+0.25*(1-np.sqrt(abs(np.sin(60*x)))))

def sinprod(x):
	return (np.sin(13*x)*np.sin(27*x)+1.)*0.5
