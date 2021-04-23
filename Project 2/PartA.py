import numpy as np

M_sol = 1.989e30        #Solar Mass
Mpc = 3.086e22          #Megaparsec Conversion
G = 6.674e-11           #Gravitational Constant Base SI


H_0 = 70 * 1000 / Mpc   #Hubble's Constant today
OmeM0 = 0.3             #Omega matter
OmeL0 = 0.7             #Omega Lambda


rho_c = (3*(H_0**2)*(Mpc**3))/(8*np.pi*G*M_sol)
                        #Critical Density

rho_m0 = OmeM0 * rho_c  #Density(matter)

rho_L = OmeL0 * rho_c   #Denstiy(Lambda)

print(f"rho_c = {rho_c} ")

z = 30                  #Redshift Start
a_1 = 1/(1+z)           #Scale factor at the z
a_0 = 1                 #Scale factor now

def timediff(a1):  
    """
    Parameters
    ----------
    a1 : Scale factor start

    Returns
    -------
    del_t : Time from when scale factor was a1 to today

    """
    del_t = (a_0 - a1)/(H_0 * a_0)
    return del_t

time_yrs = timediff(a_1) /3.1526e7 /1e9
                        #Time in millions of year since a1
                        
time_yrs0 = timediff(0) /3.1526e7 /1e9
                        #Time in billions of years since a=0 (beginning of the universe)

print(f"Time since redshift(z=30) = {time_yrs} billion years")
print(f"Time since beginning of universe(a=0)  = {time_yrs0} billion years")
print(f"Age of the Universe at (z=30) = {(time_yrs0 - time_yrs)*1e3} ")
L = 50 * Mpc
mass_50 = rho_m0 * (L**3) / M_sol

print(f"Mass in 50Mpc cube = {mass_50} Solar masses")
