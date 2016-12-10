
import numpy as np
import matplotlib.pyplot as plt
from rw import rwmodule as rw #Assumes rwmodule has been compiled with f2py to produce rw.so


def analyze_rnet(Ntime,m,X0,N0,L,Nt,display):
	"""Input variables:
	Ntime: number of time steps
    	m: number of walks
    	X0: initial node, (node with maximum degree if X0=0)
    	N0,L,Nt: recursive network parameters
    	"""
    	isample = 1
        X,Xm = rw.rwnet_omp(Ntime,m,X0,N0,L,Nt,isample,2)
        F = np.zeros([Ntime,N0+Nt])
        for t in range(int(Ntime)):
            bins = np.zeros(N0+Nt+1)
            bins[:max(X[t+1,:])+1] = np.bincount(X[t+1,:])
            F[t,:] = 1.0*bins[1:]/m
        if (display==True):
            maxloc = np.zeros(Ntime)
            maxloc = F.argmax(axis=1)+1
            plt.figure()
            plt.plot(range(1,int(Ntime+1)),maxloc,'p')
            plt.show()
        return F



def convergence_rnet(Ntime,m,X0,N0,L,Nt,display):
    	"""Input variables:
	Ntime: number of time steps
    	m: number of walks
    	X0: initial node, (node with maximum degree if X0=0)
    	N0,L,Nt: recursive network parameters
    	"""
    	isample = 1
    	M = np.zeros(m)
    	#X,Xm = rw.rwnet(Ntime,m,X0,N0,L,Nt,isample)
    	for mcount in range(1,m+1):
    	   X,Xm = rw.rwnet(Ntime,mcount,X0,N0,L,Nt,isample)
    	   M[mcount-1] = np.mean(Xm)
    	plt.figure()
    	plt.plot(range(1,m+1),M)
    	plt.show()
    	
if __name__== '__main__':
    #analyze_rnet(1e5,1e3,0,5,2,200,True)
    convergence_rnet(int(1e5),int(1e3),0,5,2,200,True)
