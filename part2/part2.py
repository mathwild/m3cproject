"""Project, part 2 by Mathilde Duverger CID:00978498"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from n1 import network as net
from p2 import flunet as fn
from time import time


def initialize(N0,L,Nt,pflag):
    """Generate network, and initial conditions
    If pflag is true, construct and return transport matrix
    """
    qmax,qnet,enet = net.generate(N0,L,Nt)
    N = (N0+Nt)
    #generate initial conditions. all nodes have S,E,C=(1,0,0)
    #except the infected node which has S,E,C=(0.1,0.05,0.05)
    init = np.zeros(3*N)
    for i in range(N):
        init[i] = 1
    #the highest degree node is infected
    infnode = qnet.argmax(axis=0)+1
    init[infnode-1] = 0.1
    init[N+infnode-1] = 0.05
    init[2*N+infnode-1] = 0.05
    if (pflag==True):
        #compute the transport matrix from the degree vector
        #and the adjacency matrix
        A = net.adjacency_matrix(N,enet)
        P = np.zeros([N,N])
        for j in range(N):
            for i in range(N):
                P[i,j] = qnet[i]*A[i,j]
            P[:,j] = P[:,j]/sum(P[:,j])
        return init,infnode,P
    return init,infnode
    
def solveFluNet(T,Ntime,a,b0,b1,g,k,w,y0,P,N,RHS=3):
    """Simulate network flu model at Ntime equispaced times
    between 0 and T (inclusive). Add additional input variables
    as needed clearly commenting what they are and how they are  
    used
    """
    #I have added the variables P the transport matrix 
    #and N the network size because they are needed
    #in the RHS.
    #I have added the variable RHS to be able to 
    #choose which RHS method we want to use when running
    #solveFluNet
    
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
    if (RHS==1):
        sol = odeint(RHSnet,y0,t,args=(a,b0,b1,g,k,w))
    if (RHS==2):
        sol = odeint(RHSnetF,y0,t,args=(a,b0,b1,g,k,w))
    if (RHS==3):
        sol = odeint(RHSnetFomp,y0,t,args=(a,b0,b1,g,k,w))
    S = sol[:,:N]
    E = sol[:,N:2*N]
    C = sol[:,2*N:3*N]
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
    init,infnode,P = initialize(N0,L,Nt,True)
    #for each w we compute Cmax,Tmax and Nmax
    for n in range((np.shape(warray))[0]):
        w = warray[n]
        #I use the Python RHS because it is the fastest so RHS=1
        t,S,E,C = solveFluNet(T,Ntime,a,b0,b1,g,k,w,init,P,N,1)
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


def performance(N0,L,Ntmax,T,Ntime,a,b0,b1,g,k,w):
    """function to analyze performance of python, fortran, and fortran+omp approaches
        Add input variables as needed, add comment clearly describing the functionality
        of this function including figures that are generated and trends shown in figures
        that you have submitted
    """
    """"I have produced two plots. The first one p24.png compares the time taken
    by solveFlunet with the three different RHS methods. We notice that the fastest
    technique is using the RHS computed in Python. This can be explained by the 
    fact that in Python we have vectorized our calculations whereas in Fortran
    we used a loop for our calculations. Secondly, we also notice that solveFlunet 
    is faster using the Fortran OMP code than using the Fortran code. This
    is expected since we use two threads with the Fortran OMP code instead of 
    only one with the regular Fortran code. The second plot p25.png provides
    the speedup of the Fortran OMP RHS compared to the Fortran RHS.
    The speedup is less than 1 for small values of Nt but increases as Nt
    increases. It goes up to 1.6 for larger sized problems, which is quite
    satisfying. We can expect that it the speedup to reach 2 for even larger
    problem sizes.
    """
    #I have added as inputs N0,L,Ntmax,T,Ntime,a,b0,b1,g,k,w to
    #be able to call initialize and solveFluNet
    Nt = Ntmax/10
    dt = Ntmax/10
    Ntu = np.zeros(10)
    
    #tpy,tf90 and tf90OMP are the times for our different RHS
    #methods
    tpy = np.zeros(10)
    tf90 = np.zeros(10)
    tf90OMP = np.zeros(10)
    
    #speedup is the speedup of the Fortran OMP RHS compared
    #to the Fortran RHS
    speedup = np.zeros(10)
    
    for i in range(10):
        #compute initial conditions
        init,infnode,P = initialize(N0,L,Nt,True)
        N = N0+Nt
        
        #time solveFluNet for the three different RHS
        t0 = time()
        solveFluNet(T,Ntime,a,b0,b1,g,k,w,init,P,N,1)
        t1 = time()
        tpy[i] = t1-t0
        t0 = time()
        solveFluNet(T,Ntime,a,b0,b1,g,k,w,init,P,N,2)
        t1 = time()
        tf90[i] = t1-t0
        t0 = time()
        solveFluNet(T,Ntime,a,b0,b1,g,k,w,init,P,N,3)
        t1 = time()
        tf90OMP[i] = t1-t0
        
        #time the two Fortran RHS methods
        t0 = time()
        fn.rhs(P,init,1,a,b0,b1,g,k,w)
        tRHS = time() - t0
        t0 = time()
        fn.rhs_omp(P,init,1,a,b0,b1,g,k,w,2)
        tRHSOMP = time() - t0
        speedup[i] = tRHS/tRHSOMP
        
        Ntu[i] = Nt
        Nt = Nt+dt
        
    plt.figure()
    plt.plot(Ntu,tpy,label='Time with Python RHS')
    plt.plot(Ntu,tf90,label='Time with Fortran RHS')
    plt.plot(Ntu,tf90OMP,label='Time with Fortran OMP RHS')
    plt.legend(loc='best')
    plt.title('performance Mathilde Duverger')
    plt.xlabel('Nt')
    plt.ylabel('Time of solveFluNet')
    
    plt.figure()
    plt.plot(Ntu,speedup,label='Speedup for RHS')
    plt.legend(loc='best')
    plt.title('performance Mathilde Duverger')
    plt.xlabel('Nt')
    plt.ylabel('Speedup of the Fortran OMP code compared to the Fortran code')
    plt.show()
    


if __name__ == '__main__':            
   a,b0,b1,g,k,w = 45.6,750.0,0.5,73.0,1.0,0.1
   warray = np.array([0,1e-2,1e-1,0.2,0.5,1.0])
   Cmax,Tmax,Nmax = analyze(5,2,500,2,100,a,b0,b1,g,k,0.1,warray,True)
   performance(5,2,500,2,100,a,b0,b1,g,k,w)

