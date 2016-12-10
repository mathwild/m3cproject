import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from p3 import sync #assumes that fortran module sync has been compiled with f2py as p3.so



def oscillator(Nt,T,N,c,mu=1,s):
    """Simulate fully-coupled network of oscillators
    Compute simulation for Nt time steps between 0 and T (inclusive)
    input: N: number of oscillators
           c: coupling coefficient
           mu,sigma: distribution parameters for omega_i
    """
    
    
    

if __name__ == '__main__':
    n,c,m,s = 101,10.0,1.0,0.1
    Nt,T = 500,100
    t,omega,theta0,theta,order = kserial(Nt,T,n,c,m,s)
