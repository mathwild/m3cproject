import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from p3 import sync #assumes that fortran module sync has been compiled with f2py as p3.so



def oscillator(Nt,T,N,c,mu,s):
    """Simulate fully-coupled network of oscillators
    Compute simulation for Nt time steps between 0 and T (inclusive)
    input: N: number of oscillators
           c: coupling coefficient
           mu,sigma: distribution parameters for omega_i
    """
    sync.ntotal = N
    sync.c = c
    sync.w = np.random.normal(mu,s,N)
    theta0  = np.random.uniform(0,2*np.pi,N)
    dt = T/float(Nt)
    theta,order = sync.rk4(0,theta0,dt,Nt)
    return theta,order
    
    

if __name__ == '__main__':
    n,c,m,s = 101,10.0,1.0,0.1
    Nt,T = 500,100
    #theta,order = oscillator(Nt,T,n,c,m,s)
    #plt.figure()
    #plt.plot(order)
    #plt.figure()
    #plt.plot(np.mod(np.transpose(theta),2*np.pi))
    #plt.show()
    plt.figure()
    plt.scatter(range(101),np.mod(np.loadtxt('theta.dat'),2*np.pi))
    plt.figure()
    plt.plot(range(500),np.loadtxt('order.dat'))
    plt.show()