import math
import numpy as np
from constants import *

#Secondary Functions included in the ODE's

def Rho(mu,P,T):
    """
    Density with respect to the radius

    Parameters
    ----------
    mu : Mean Molecular Weight
    
    P : Pressure
    
    T : Temperature

    Returns
    -------
    
    Density (Dimensionless)

    """
    rho = ((P_sol*P) - (a*((T_sol*T)**4)/3))*(mu*Mh)/(k*T_sol*T)
    
    return rho/rho_sol
    
def Eps(mu,X,Y,Z,P,T):
    """
    Total Energy released per kg per s by all nuclear reactions

    Parameters
    ----------
    X : Hydrogen composition
    
    Y : Helium composition
    
    Z : Metallic composition
    

    Returns
    -------
    Epsilon (Dimensionless)

    """
    rho = Rho(mu, P, T)*rho_sol
    T6 = (T_sol*T)/1e6
    
    psi_pp = 1 + 1.412e8*(1/X - 1)*np.exp(-49.98*(T6**(-1/3)))
    Cpp = 1 + 0.012*(T6**(1/3)) + 0.0109*(T6**(2/3)) + 0.000938*T6
    
    eps_pp = 0.241 * rho * X * X * fpp * psi_pp * Cpp * (T6**(-2/3)) * np.exp(-33.80*(T6**(-1/3)))
    
    
    Xcno = Z/2
    
    Ccno = 1 + 0.0027*(T6**(1/3)) - 0.00778*(T6**(2/3)) - 0.000149*T6
    
    eps_cno = 8.67e20 * rho * X * Xcno * Ccno * (T6**(-2/3)) * np.exp(-152.28*(T6**(-1/3)))
    
    
    T8 = (T_sol*T)/1e8
    
    eps_3a = 50.9 * rho * rho * Y * Y * Y * (T8**(-3)) * f3a * np.exp(-44.027*(T8**(-1)))
    
    eps = eps_pp + eps_cno + eps_3a
    
    return eps/eps_sol

def k_c(mu, X, Y, Z, P, T):
    """
    Opacities

    Parameters
    ----------
    As listed above.

    Returns
    -------
    kappa

    """
    rho = Rho(mu, P, T)*rho_sol
    gbf_t = 1 / (0.708 * ((rho * (1+X))**(1/5)))
    
    kcbf = 4.34e21 * (gbf_t) * Z * (1 + X) * rho /((T_sol*T)**3.5)
    
    
    kcff = 3.68e18 * gff * (1 - Z) * (1 + X) * rho / ((T_sol*T)**3.5)
    
    kces = 0.02 * (1 + X)
    
    if 3000 <= T_sol*T <= 6000 and 1e-7 <= rho*rho_sol <= 1e-2 and 0.001 < Z < 0.03:
        kH = 7.9e-34 * (Z/0.02) * (rho**(1/2)) * ((T_sol*T)**9)
    else:
        kH = 0
        
    kc = kcbf + kcff + kces + kH
    
    return kc/kc_sol


