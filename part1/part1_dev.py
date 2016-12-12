
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



def convergence_rnet(Ntime,m,X0,N0,L,Nt):
    	"""Input variables:
	Ntime: number of time steps
    	m: number of walks
    	X0: initial node, (node with maximum degree if X0=0)
    	N0,L,Nt: recursive network parameters
    	"""
    	isample = 1
    	X,Xm = rw.rwnet(Ntime,m,X0,N0,L,Nt,isample)
    	Xmcount = np.zeros(m)
    	Xmcount[m-1] = Xm[Ntime-1]
    	marray= np.zeros(m)
    	marray[m-1] = m
    	for mcount in range(1,m):
    	   Xt = X[:,:mcount]
    	   bins = np.count_nonzero(Xt[Ntime,:]==Xt[0,0])
    	   Xmcount[mcount-1] = 1.0*bins/mcount
    	   marray[mcount-1] = mcount
    	plt.figure()
    	plt.plot(marray,Xmcount)
    	plt.show()
    	
if __name__== '__main__':
    #analyze_rnet(1e5,1e3,0,5,2,200,True)
    convergence_rnet(1e5,1000,0,5,2,200)
