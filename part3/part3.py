"""Project part 3 by Mathilde Duverger CID: 00978498
"""

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
    t = np.linspace(0,T,Nt)
    sync.ntotal = N
    sync.c = c
    omega = np.random.normal(mu,s,N)
    sync.w = omega
    theta0  = np.random.uniform(0,2*np.pi,N)
    dt = T/float(Nt)
    theta,order = sync.rk4(0,theta0,dt,Nt)
    return t,omega,theta0,theta,order
    

if __name__ == '__main__':
    n,c,m,s = 101,10.0,1.0,0.1
    Nt,T = 500,100
    t,omega,theta0,theta,order = oscillator(Nt,T,n,c,m,s)