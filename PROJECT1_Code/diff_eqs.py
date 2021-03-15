import math
import numpy as np 
import bterms as b
import mod_rk4 as rk4
from constants import *
import matplotlib.pyplot as plt

#ODE's


def funcP(mu, X, Y, Z, r, P, M, L ,T):
    """
    Pressure ODE dP/dR - Dimensionless
    
    Parameters
    ----------
    mu : Mean Molecular Weight
    
    mu : Mean Molecular Weight
    
    X : Hydrogen composition
    
    Y : Helium Composition
    
    Z : Metallic composition
    
    r : Radius
    
    P : Pressure 
    
    M : Mass
    
    L : Luminosity
    
    T : Temperature

    Returns
    -------
    Pressure

    """
    rho = b.Rho(mu, P, T)
    
    p = c1 * M * rho /(r*r)
    
    return p

def funcM(mu, X, Y, Z, r, P, M, L ,T):
    """
    Mass ODE dMr/dr - Dimensionless

    Parameters
    ----------
    Same as above

    Returns
    -------
    Mass

    """
    rho = b.Rho(mu, P, T)
    
    m = c2 * rho * r * r
    return m

def funcL(mu, X, Y, Z, r, P, M, L ,T):
    """
    Luminosity ODE dL/dr - Dimensionless

    Parameters
    ----------
    rho - Density
    eps - Epsilon : Function of Energy released from Nuclear reaction

    Returns
    -------
    luminosity

    """
    rho =  b.Rho(mu, P, T)
    eps = b.Eps(mu, X, Y, Z, P, T)
    
    l = c3 * rho * eps * r * r
    return l

def funcT(mu, X, Y, Z, r, P, M, L ,T):
    """
    Radiative Temperature ODE - Dimensionless

    Parameters
    ----------
    kc - Opacity
    
    Returns
    -------
    Temperature

    """
    rho = b.Rho(mu, P, T)
    kc = b.k_c(mu, X, Y, Z, P ,T)
    
    t = c4 * kc * rho * L /((T**3)*(r**2))
    
    return t

def funcT2(mu, X, Y, Z, r, P, M, L ,T):
    """
    Convective Temperature ODE - Dimensionless

    Parameters
    ----------
    As above

    Returns
    -------
    Temperature

    """
    
    t2 = c5 * (M * mu)/ (r * r)
    return t2


def PMLT_plot():
    """
    Plots P, M, L or T against radius

    """
    Xtest = 0.711
    Ytest = 0.254
    Ztest = 1 - Xtest - Ytest
    htest = 0.0001
    mutest = (2*Xtest + (3*Ytest/4) + Ztest/2)**(-1)
    rtest, Ptest, Mtest, Ltest, Ttest = rk4.solveODE(30, mutest, Xtest, Ytest, Ztest, htest, funcP, funcM, funcL, funcT, funcT2)
    
    plt.plot(rtest,Ltest)
    plt.show()

#PMLT_plot()


def multistars_comp(N):
    """
    Models N number of stars

    Parameters
    ----------
    X - Array of 20 random hydrogen composition values between 0.7 and 0.87
    
    Y - Array of 20 random hydrogen composition values between 0.1 and 0.26
    
    Z - Metallic composition

    Returns
    -------
    None.

    """
    
    X = np.arange(0.7,0.87, (0.87-0.7)/N)
    Y_inv = np.arange(0.1,0.26, (0.26-0.1)/N)
    Y = np.flip(Y_inv)
    Z = 1 - X - Y
    
    mstop = 30
    
    h = 0.00001
    
    r = np.zeros(N)
    P = np.zeros(N)
    M = np.zeros(N)
    L = np.zeros(N)
    T = np.zeros(N)
    
    #Takes runs the ODE N times, takes the last values of each ODE
    for i in range(N):
        mu = (2*X[i] + (3*Y[i]/4) + Z[i]/2)**-1
        r_end,P_end,M_end,L_end,T_end = rk4.solveODE(mstop,mu, X[i], Y[i], Z[i], h, funcP, funcM, funcL, funcT, funcT2)
        
        r[i] = r_end[-1]
        P[i] = P_end[-1]
        M[i] = M_end[-1]
        L[i] = L_end[-1]
        T[i] = L[i]/(4*np.pi*r[i]*r[i]*sigma*T_sol )
    
    print(F"LENGTH = {len(r)}")
    values = np.column_stack((r,P,M,L,T))
    np.savetxt("comptest.txt",values)       #Saves data to a txt file
    
    return

multistars_comp(20)


