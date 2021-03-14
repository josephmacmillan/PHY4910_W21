# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 16:21:12 2021

@author: Darryen
"""
import numpy as np
import mpmath as mp


# Constants
G = 6.67E-11
c = 3E8
h = 6.63E-34
kB = 1.38E-23
mE = 9.11E-31
mH = 1.67E-27
mSun = 2E30
LSun = 3.8E26
rSun = 7E8
TSun = 1E7
a = 7.5E-16

def mu(X,Y,Z):
    return 1/(2*X + 0.75 * Y + 0.5*Z)


def kbf(X,Z,rho, T):
    '''
    Parameters
    ----------
    X : Float
        Hydrogen mass fraction
    Z : Float
        Metal mass fraction
    rho : Float
        The density
    T : Float
        Temperature

    Returns
    -------
    kbf : Float
        The opacity from the bound-free processes
    '''
    if np.isnan((4.34E21)*(0.708*(rho*(1+X))**(1/5))**(-1)*Z*(1+X)*rho/(T ** 3.5)) != True:
        print("KBF:", (4.34E21)*(0.708*(rho*(1+X))**(1/5))**(-1)*Z*(1+X)*rho/(T ** 3.5))
    return (4.34E21)*(0.708*(rho*(1+X))**(1/5))**(-1)*Z*(1+X)*rho/(T ** 3.5)


def kff(X,Z,rho, T):
    '''
    Parameters
    ----------
    X : Float
        Hydrogen mass fraction
    Z : Float
        Metal mass fraction
    rho : Float
        The density
    T : Float
        Temperature

    Returns
    -------
    kff : Float
        The opacity from the free-free processes
    '''
    if np.isnan(((3.68E18)*(1-Z)*(1+X)*rho)/(T ** 3.5)) != True:
        print("KFF: ", ((3.68E18)*(1-Z)*(1+X)*rho)/(T ** 3.5))
    return ((3.68E18)*(1-Z)*(1+X)*rho)/(T ** 3.5)

def kes(X):
    '''
    X : Float
        Hydrogen mass fraction
    Returns
    -------
    kes : Float
        The opacity from the electron-scattering processes
    '''
    return (0.02)*(1+X)

def kh(Z, rho, T):
    '''
    Parameters
    ----------
    Z : Float
        Metal mass fraction
    rho : Float
        The density
    T : Float
        Temperature

    Returns
    -------
    kh : Float
        The opacity from the H-ion processes
    '''
    if np.isnan((7.94E-34)*(Z/0.02)*(rho ** (1/2))*T**9) != True:
        print("KH: ", (7.94E-34)*(Z/0.02)*(rho ** (1/2))*T**9)
    if T < 6000 and rho < 1E-2 and Z < 0.03:
        return (7.94E-34)*(Z/0.02)*(rho ** (1/2))*T**9
    else:
        return 0


def k(X, Z, rho, T):
    '''
    Parameters
    ----------
    kbf : Function
        Bound-free processes
    kff : Function
        Free-free processes
    kes : Function
        Electron-scattering processes
    kh : Function
        H- ion processes

    Returns
    -------
    k : Float
        The opacity of the star

    '''
    #print('Opacity: ', kbf(X,Z,rho, T) + kff(X,Z,rho, T) + kes(X) + kh(Z, rho, T))
    return kbf(X,Z,rho, T) + kff(X,Z,rho, T) + kes(X) + kh(Z, rho, T)


def epp(X, rho, T):
    '''
    Parameters
    ----------
    X : Float
        Hydrogen mass fraction
    rho: Float
        Density
    T : Float
        Temperature
    Returns
    -------
    epp : Float
        Energy per kg per s by proton-proton chain
    '''
    T6 = T / 1E6
    psipp = 1 + 1.412E8*(1/X - 1)*np.exp(-49.98*(T6 ** (-1/3)))
    cpp = 1 + 0.012*(T6 ** (1/3)) + 0.0109*(T6 ** (2/3)) + 0.000938*T6
    
    if np.isnan((0.241)*rho*X**2 *psipp*cpp*(T6 ** (-2/3))*np.exp(-33.80*(T6 ** (-1/3)))) != True:
    
        print("EPP: ", (0.241)*rho*X**2 *psipp*cpp*(T6 ** (-2/3))*np.exp(-33.80*(T6 ** (-1/3))))
    
    return (0.241)*rho*X**2 *psipp*cpp*(T6 ** (-2/3))*np.exp(-33.80*(T6 ** (-1/3)))


def ecno(X,Z,rho, T):
    '''
    Parameters
    ----------
    X : Float
        Hydrogen mass fraction
    Z : Float
        Metal mass fraction
    rho: Float
        Density
    T : Float
        Temperature
    Returns
    -------
    ecno : Float
        Energy per kg per s by CNO cycle
    '''
    T6 = T/ 1E6
    xcno = Z/2
    ccno = 1 + 0.0037*(T6 ** (1/3)) - 0.00778*(T6 ** (2/3)) - 0.000149*T6
    
    if np.isnan((8.67E20)*rho*ccno*xcno*X*(T6 ** (-2/3))*np.exp(-152.28*(T6 ** (-1/3)))) != True:
        print("ECNO :", (8.67E20)*rho*ccno*xcno*X*(T6 ** (-2/3))*np.exp(-152.28*(T6 ** (-1/3))))
    return (8.67E20)*rho*ccno*xcno*X*(T6 ** (-2/3))*np.exp(-152.28*(T6 ** (-1/3)))

def e3a(Y,rho,T):
    '''
    Parameters
    ----------
    Y : Float
        Helium mass fraction
    rho : Float
        Density
    T : Float
        Temperature

    Returns
    -------
    e3a : Float
        Energy per kg per s by 3 - alpha processes
    '''
    T8 = T / 1E8
    if np.isnan((50.9)*rho**2 * Y**3 * (T8 ** (-3)) * np.exp(-44.027*T8**(-1))) != True :
        print("E3A: ", (50.9)*rho**2 * Y**3 * (T8 ** (-3)) * np.exp(-44.027*T8**(-1)))
    return (50.9)*rho**2 * Y**3 * (T8 ** (-3)) * np.exp(-44.027*T8**(-1))

def e(X, Y, Z, rho, T):
    '''
    Parameters
    ----------
    epp : Function
        Energy per kg per s by proton-proton chain
    ecno : Function
        Energy per kg per s by CNO cycle
    e3a : Function
        Energy per kg per s by 3-alpha processes

    Returns
    -------
    e : Float
        The total energy per kg per s by all nuclear reactions
    '''
    #print("Energy: ", epp(X, rho, T) + ecno(X,Z,rho, T) + e3a(Y,rho,T))
    return epp(X, rho, T) + ecno(X,Z,rho, T) + e3a(Y,rho,T)


def density(X,Y,Z, p, T):
    '''
    Parameters
    ----------
    p : Float
        Pressure of star
    T : Float
        Temperature
    X : Float
        Hydrogen mass fraction
    Y : Float
        Helium mass fraction
    Z : Float
        Metal mass fraction

    Returns
    -------
    rho : Float
        Density 
    '''
    print('Density: ' , (p - (1/3)*a*T**3)*((mu(X,Y,Z)*mH)/(kB)))
    return (p - (1/3)*a*T**3)*((mu(X,Y,Z)*mH)/(kB))