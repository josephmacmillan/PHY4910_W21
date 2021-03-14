# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 20:53:55 2021

@author: Darryen
"""
import numpy as np
import RK4 as RK
import HelperFunctions as HF


# Mass fractions from website
# http://burro.case.edu/Academics/Astr221/StarPhys/stellarint.html
X = 0.7
Y = 0.28
Z = 1 - X - Y

def dp(r, p, M, L, T):
    '''
    Parameters
    ----------
    r : Float
        Radius
    p : Float
        Pressure
    M : Float
        Mass
    L : Float
        Luminosity 
    T : Float
        Temperature

    Returns
    -------
    dp/dr : Float
            The change in pressure wrt radius.
    '''
    rho = HF.density(X,Y,Z, p,T)
    print('dp:',(- HF.G * M * rho)/(r**2))
    return (- HF.G * M * rho)/(r**2)

def dM(r, p, M, L, T):
    '''
    Parameters
    ----------
    r : Float
        Radius
    p : Float
        Pressure
    M : Float
        Mass
    L : Float
        Luminosity 
    T : Float
        Temperature

    Returns
    -------
    dM/dr : Float
            The change in mass wrt radius.
    '''
    rho = HF.density(X,Y,Z, p,T)
    print('dM: ',4*np.pi*r**2 * rho)
    return 4*np.pi*r**2 * rho

def dL(r,p,M,L,T):
    '''
    Parameters
    ----------
    r : Float
        Radius
    p : Float
        Pressure
    M : Float
        Mass
    L : Float
        Luminosity 
    T : Float
        Temperature

    Returns
    -------
    dL/dr : Float
            The change in luminosity wrt radius.
    '''
    rho = HF.density(X,Y,Z, p,T)
    e = HF.e(X, Y, Z, rho, T)
    print('dL:',(4*np.pi * r **2 * rho * e))
    return (4*np.pi * r **2 * rho * e)

def dT(r,p,M, L, T):
    '''
    Parameters
    ----------
    r : Float
        Radius
    p : Float
        Pressure
    M : Float
        Mass
    L : Float
        Luminosity 
    T : Float
        Temperature

    Returns
    -------
    dT/dr : Float
            The change in temperature wrt radius.
    '''
    rho = HF.density(X,Y,Z, p,T)
    k = HF.k(X, Z, rho, T)
    print('dT:',(-3*k*rho * L)/(4*HF.a*HF.c*T**3 * 4 * np.pi * r ** 2))
    return (-3*k*rho * L)/(16*HF.a*HF.c*T**3  * np.pi * r ** 2)

    

def SymStatStar(end, stepsize, p0, M0, L0, T0):
    '''
    Parameters
    ----------
    M : Float
        Total mass of star
    X : Float
        Hydrogen mass fraction 
    Y : Float
        Helium mass fraction
    Z : Float
        Metal mass fraction
    
    Returns
    -------
    eta : Float
        The surface radius of the star
    T : Float
        The temperature of the star
    L : Float
        The luminosity of the star
    rho : Float
        The density of the star
    p : Float
        The pressure of the star
    '''
    r, p, M, L, T = RK.ode_rk4(0.1*stepsize, end, stepsize, p0, M0, L0, T0, dp, dM, dL, dT)
    rho = HF.density(X,Y,Z, p,T)
    
    print(f"Radius: {r[-1]} \n Temperature: {T[-1]} \n Luminosity: {L[-1]} \n Density: {rho[-1]} \n Pressure: {p[-1]}")
    
    return r, T, L, rho, p, M

radius, T, L, rho, p, M = SymStatStar(2500, 100, 2.5E16, 0,0,1.5E7)
    
