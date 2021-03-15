import numpy as np
from constants import*
import math
import bterms as b


def solveODE(mstop, mu, X, Y, Z, h, func_E, func_F, func_G, func_H, func_I): 
    """
    Solves 4 ODe''

    Parameters
    ----------
    mstop : Stops if the solver approaches a specific mass
    
    mu : Mean Molecular Weight
    
    X : Hydrogen composition
    
    Y : Helium Composition
    
    Z : Metallic composition
    
    h : steps for solver
    
    func_E : Pressure ODE
    
    func_F : Mass ODE
    
    func_G : Luminosity ODE
    
    func_H : Radiative Transport ODE
    
    func_I : Convective Transport ODE

    Returns
    -------
    r, P, M, L, T arrays
    """
    
    #initialize the arrays for each variable
    r = np.arange(1e-5,2,h)
    N = len(r)
    P = np.zeros(N)
    M = np.zeros(N)
    L = np.zeros(N)
    T = np.zeros(N)
    
    #set initial conditions
    P[0] = 2.54
    M[0] = 0
    L[0] = 0
    T[0] = 1.52
    
    #Condition for choosing between radiative and convective transport
    condition = True
    
    
    # default to runge-kutta
    for i in range(0,N-1):
        
        # Loops to find new P, M, L,and T based on functions and ji, ki, li, mi values
        j1 = h * func_E(mu, X, Y, Z, r[i],P[i],M[i],L[i],T[i])
        k1 = h * func_F(mu, X, Y, Z, r[i],P[i],M[i],L[i],T[i])
        l1 = h * func_G(mu, X, Y, Z, r[i],P[i],M[i],L[i],T[i])
        if condition == True:
            m1 = h * func_H(mu, X, Y, Z, r[i],P[i],M[i],L[i],T[i])
        else:
            m1 = h * func_I(mu, X, Y, Z, r[i],P[i],M[i],L[i],T[i])
        
       
        
            
        j2 = h * func_E(mu, X, Y, Z, r[i] + h/2, P[i] + j1/2, M[i] + k1/2, L[i] + l1/2, T[i] + m1/2)
        k2 = h * func_F(mu, X, Y, Z, r[i] + h/2, P[i] + j1/2, M[i] + k1/2, L[i] + l1/2, T[i] + m1/2)
        l2 = h * func_G(mu, X, Y, Z, r[i] + h/2, P[i] + j1/2, M[i] + k1/2, L[i] + l1/2, T[i] + m1/2)
        if condition == True:
            m2 = h * func_H(mu, X, Y, Z, r[i] + h/2, P[i] + j1/2, M[i] + k1/2, L[i] + l1/2, T[i] + m1/2)
        else:
            m2 = h * func_I(mu, X, Y, Z, r[i] + h/2, P[i] + j1/2, M[i] + k1/2, L[i] + l1/2, T[i] + m1/2)
            
      
            
        j3 = h * func_E(mu, X, Y, Z, r[i] + h/2, P[i] + j2/2, M[i] + k2/2, L[i] + l2/2, T[i] + m2/2)
        k3 = h * func_F(mu, X, Y, Z, r[i] + h/2, P[i] + j2/2, M[i] + k2/2, L[i] + l2/2, T[i] + m2/2)
        l3 = h * func_G(mu, X, Y, Z, r[i] + h/2, P[i] + j2/2, M[i] + k2/2, L[i] + l2/2, T[i] + m2/2)
        if condition == True:
            m3 = h * func_H(mu, X, Y, Z, r[i] + h/2, P[i] + j2/2, M[i] + k2/2, L[i] + l2/2, T[i] + m2/2)
        else:
            m3 = h * func_I(mu, X, Y, Z, r[i] + h/2, P[i] + j2/2, M[i] + k2/2, L[i] + l2/2, T[i] + m2/2)
            
        
            
            
        j4 = h * func_E(mu, X, Y, Z, r[i] + h, P[i] + j3, M[i] + k3, L[i] + l3, T[i] + m3)
        k4 = h * func_F(mu, X, Y, Z, r[i] + h, P[i] + j3, M[i] + k3, L[i] + l3, T[i] + m3)
        l4 = h * func_G(mu, X, Y, Z, r[i] + h, P[i] + j3, M[i] + k3, L[i] + l3, T[i] + m3)
        if condition == True:
            m4 = h * func_H(mu, X, Y, Z, r[i] + h, P[i] + j3, M[i] + k3, L[i] + l3, T[i] + m3)
        else:
            m4 = h * func_I(mu, X, Y, Z, r[i] + h, P[i] + j3, M[i] + k3, L[i] + l3, T[i] + m3)
               
        
            
        #use values (derivatives) to compute Pi+1 etc (following points)
        P[i+1] = P[i] + j1/6 + j2/3 + j3/3 + j4/6
        M[i+1] = M[i] + k1/6 + k2/3 + k3/3 + k4/6
        L[i+1] = L[i] + l1/6 + l2/3 + l3/3 + l4/6
        T[i+1] = T[i] + m1/6 + m2/3 + m3/3 + m4/6
        
        #checks the temperature gradtient of adjacent points for convective take over
        rad =(np.log(P[i+1] - np.log(P[i])))/(np.log(T[i+1]) - np.log(T[i]))
        #(T[i]*func_E(mu, X, Y, Z, r[i],P[i],M[i],L[i],T[i]))/(P[i]*func_H(mu, X, Y, Z, r[i],P[i],M[i],L[i],T[i]))#np.abs(func_H(mu, X, Y, Z, r[i],P[i],M[i],L[i],T[i]))
        con = gamma/(gamma -1)
        if rad > con:
            condition = False
        else:
            condition = True

        #returns r, P,M,L,T when temperature, and pressure goes to zero, or solver returns nan
        if T[i+1] <= 1e-15 or P[i+1]<=1e-10 or math.isnan(T[i+1])==True: #M[i+1] >= mstop and
            return r[0:i], P[0:i], M[0:i], L[0:i], T[0:i]
    



