import numpy as np

H0 = 68*(1000/3.086e+22)
omega_M = 0.3089
tot_density = 1 
omega_L = tot_density - omega_M
G = 6.673e-11 

#calculate critcal density today (rho_c,0). Units are in kg/m^3 
rho_c0 = (3*(H0**2))/(8*np.pi*G)
print('The current critical density of the universe is', rho_c0,'kg/m^3')

#calculate matter density today rho_m,0
rho_m0 = rho_c0 * omega_M
print('The current matter density is', rho_m0,'kg/m^3')


# Calculating look back times

z = 30
a = 1/(1+z)


def lookback_time(a1):
    """
    
    Parameters
    ----------
    a1 : int 
        the initial scale factor.

    Returns
    -------
    t : the look back time. 

    """
    a0 = 1
    
    t = (a0 - a1)/(H0*a0) 
    return t 

age_of_univ = (lookback_time(0))/(3.154e+7*1e+9)
print('The current age of the universe is',age_of_univ,'billion years')

age_start_sim = (lookback_time(a))/(3.154e+7*1e+6)
print('The age of the universe at the start of the simulation is',age_start_sim,'million years')

#Calculate the total mass of matter in a cube with side lenghts of 50 Mpc. 
L = 50 #Mpc 

mass = (L**3)*(rho_c0*(3.086e+22**3))
print('The total mass inside the box is',mass ,'kg')