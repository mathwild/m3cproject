"""Project, part 2"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from n1 import network as net
from p2 import flunet as fn
#from time import clock,time


def initialize(N0,L,Nt,pflag):
    """Generate network, and initial conditions
    If pflag is true, construct and return transport matrix
    """
    qmax,qnet,enet = net.generate(N0,L,Nt)
    N = (N0+Nt)
    init = np.zeros(3*N)
    for i in range(N):
        init[i] = 1
    infnode = qnet.argmax(axis=0)+1
    init[infnode-1] = 0.1
    init[N+infnode-1] = 0.05
    init[2*N+infnode-1] = 0.05
    if (pflag==True):
        A = net.adjacency_matrix(N,enet)
        P = np.zeros([N,N])
        for j in range(N):
            for i in range(N):
                P[i,j] = qnet[i]*A[i,j]
            P[:,j] = P[:,j]/sum(P[:,j])
        return init,infnode,P
    return init,infnode
    
def solveFluNet(N0,L,Nt,T,Ntime,a,b0,b1,g,k,w):
    """Simulate network flu model at Ntime equispaced times
    between 0 and T (inclusive). Add additional input variables
    as needed clearly commenting what they are and how they are  
    used
    """
    init,infnode,P = initialize(N0,L,Nt,True)
    N = N0+Nt
    #add input variables to RHS functions if needed
    def RHSnet(y,t,a,b0,b1,g,k,w):
        """RHS used by odeint to solve Flu model"""
        S = y[:N]
        E = y[N:2*N]
        C = y[2*N:3*N]
        b = b0 + b1*(1+np.cos(2*np.pi*t))
        dy = np.zeros(3*N)
        dy[:N]= k*(1-S)-b*C*S+w*np.dot(P,S)-w*S
        dy[N:2*N]= b*C*S-(k+a)*E+w*np.dot(P,E)-w*E
        dy[2*N:3*N]= a*E-(g+k)*C+w*np.dot(P,C)-w*C
        return dy
    
    def RHSnetF(y,t,a,b0,b1,g,k,w):
        """RHS used by odeint to solve Flu model"
        Calculations carried out by fn.rhs
        """
        dy = fn.rhs(P,y,t,a,b0,b1,g,k,w)
        return dy
        
    def RHSnetFomp(y,t,a,b0,b1,g,k,w):
        """RHS used by odeint to solve Flu model
        Calculations carried out by fn.rhs_omp
        """
        dy = fn.rhs_omp(P,y,t,a,b0,b1,g,k,w,2)
        return dy

    #Add code here and to RHS functions above to simulate network flu model
    t = np.linspace(0,T,Ntime)
    sol = odeint(RHSnet,init,t,args=(a,b0,b1,g,k,w))
    S = sol[:,:N]
    E = sol[:,N:2*N]
    C = sol[:,2*N:3*N]
    #dy1 = RHSnet(init,1,a,b0,b1,g,k,w)
    #dy2 = RHSnetF(init,1,a,b0,b1,g,k,w)
    #print 'shapedy1', np.shape(dy1)
    #print 'shapedy2', np.shape(dy2)
    #print 'dy1', dy1
    #print 'dy2', dy2
    return t,S,E,C


def analyze(N0,L,Nt,T,Ntime,a,b0,b1,g,k,threshold,warray,display=False):
    """analyze influence of omega on: 
    1. maximum contagious fraction across all nodes, Cmax
    2. time to maximum, Tmax
    3. maximum fraction of nodes with C > threshold, Nmax    
    Input: N0,L,Nt: recursive network parameters 
           T,Ntime,a,b0,b1,g,k: input for solveFluNet
           threshold: use to compute  Nmax
           warray: contains values of omega
    """
    Cmax = np.zeros((np.shape(warray))[0])
    Tmax = np.zeros((np.shape(warray))[0])
    Nmax = np.zeros((np.shape(warray))[0])
    N = (N0+Nt)
    for n in range((np.shape(warray))[0]):
        w = warray[n]
        t,S,E,C = solveFluNet(N0,L,Nt,T,Ntime,a,b0,b1,g,k,w)
        Cmax[n] = np.amax(C)
        u = [(index) for index,row in enumerate(C) if np.amax(C) in row]
        Tmax[n] = u[0] + 1
        NC = np.zeros(Ntime)
        for i in range(Ntime):
            count = 0
            for j in range(N):
                if (C[i,j]>threshold):
                    count = count+1
            NC[i] = 1.0*count/N   
        Nmax[n] = max(NC)
    if (display==True):
        plt.figure()
        plt.plot(warray,Cmax,label='Cmax')
        plt.xlabel('w')
        plt.legend(loc='best')
        plt.title('analyze Mathilde Duverger')
        plt.figure()
        plt.plot(warray,Tmax,label='Tmax')
        plt.xlabel('w')
        plt.legend(loc='best')
        plt.title('analyze Mathilde Duverger')
        plt.figure()
        plt.plot(warray,Nmax,label='Nmax')
        plt.xlabel('w')
        plt.legend(loc='best')
        plt.title('analyze Mathilde Duverger')
        plt.show()
    return Cmax,Tmax,Nmax
    

def visualize(enet,C,threshold):
    """Optional, not for assessment: Create figure showing nodes with C > threshold.
    Use crude network visualization algorithm from homework 3
    to display nodes. Contagious nodes should be red circles, all
    others should be black"""
    return None


def performance():
    """function to analyze performance of python, fortran, and fortran+omp approaches
        Add input variables as needed, add comment clearly describing the functionality
        of this function including figures that are generated and trends shown in figures
        that you have submitted
    """


if __name__ == '__main__':            
   a,b0,b1,g,k,w = 45.6,750.0,0.5,73.0,1.0,0.1
   #t,S,E,C = solveFluNet(5,2,2,10,100,a,b0,b1,g,k,w)
   #print 'S=',S
   #print 'E=', E
   #print 'C=', C
   warray = np.array([0,1e-2,1e-1,0.2,0.5,1.0])
   #warray = np.array([0.1])
   Cmax,Tmax,Nmax = analyze(5,2,500,2,100,a,b0,b1,g,k,0.1,warray,True)
   #print 'Cmax=', Cmax
   #print 'Tmax=', Tmax
   #print 'Nmax=', Nmax


