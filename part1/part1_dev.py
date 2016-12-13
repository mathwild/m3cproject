
import numpy as np
import matplotlib.pyplot as plt
from rw import rwmodule as rw #Assumes rwmodule has been compiled with f2py to produce rw.so


def analyze_rnet(Ntime,m,X0,N0,L,Nt,display=False):
	"""Input variables:
	Ntime: number of time steps
    	m: number of walks
    	X0: initial node, (node with maximum degree if X0=0)
    	N0,L,Nt: recursive network parameters
    	"""
    	isample = 1
        X,Xm = rw.rwnet_omp(Ntime,m,X0,N0,L,Nt,isample,2)
        F = np.zeros([Ntime,N0+Nt])
        #compute F the fraction of the m walkers
        #at node n at time t
        for t in range(int(Ntime)):
            bins = np.zeros(N0+Nt+1)
            bins[:max(X[t+1,:])+1] = np.bincount(X[t+1,:])
            F[t,:] = 1.0*bins[1:]/m
        if (display==True):
            maxloc = np.zeros(Ntime)
            maxloc = F.argmax(axis=1)+1
            #plot the node with the greatest number of walkers
            #at each step
            plt.figure()
            plt.plot(range(1,int(Ntime+1)),maxloc,'p')
            plt.title('analyze_rnet Mathilde Duverger')
            plt.xlabel('Step')
            plt.ylabel('Node with greatest walkers')
            plt.show()
        return F



def convergence_rnet(Ntime,m,X0,N0,L,Nt,display=False):
        """This code computes the fraction of walkers at the initial
        node in the final timestep for different number of walkers, m.
        From the plot, we notice that the fraction of walkers firstly increases
        converges to a value between 0.04 and 0.05. """
    	"""Input variables:
	Ntime: number of time steps
    	m: number of walks
    	X0: initial node, (node with maximum degree if X0=0)
    	N0,L,Nt: recursive network parameters
    	"""
    	isample = 1
    	X,Xm = rw.rwnet(Ntime,m,X0,N0,L,Nt,isample)
    	F = np.zeros(m)
    	F[m-1] = Xm[Ntime-1]
    	for mcount in range(1,m):
    	   Xt = X[:,:mcount]
    	   bins = np.count_nonzero(Xt[Ntime,:]==Xt[0,0])
    	   F[mcount-1] = 1.0*bins/mcount
    	if (display==True):
    	   plt.figure()
    	   plt.plot(range(1,m+1),F,label='Fraction of walkers at the initial node in the final timestep')
    	   plt.title('convergence_rnet Mathilde Duverger')
    	   plt.xlabel('m')
    	   plt.legend(loc='best')
    	   plt.show()
    	
if __name__== '__main__':
    analyze_rnet(1e5,1e3,0,5,2,200,True)
    #convergence_rnet(1e5,1000,0,5,2,200,True)
