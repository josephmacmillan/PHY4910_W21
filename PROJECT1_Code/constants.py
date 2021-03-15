import numpy as np
import math


G = 6.67430e-11         #Gravitational Constant
c = 3e8                 #Speed of light
k = 1.38064852e-23      #Boltzmann constant
gamma = 5/3             #Ideal Gas capacity Ratio
Mh = 1.67e-27           #Mass of hydrogen
sigma = 5.67037e-8      #Stephen-Boltzmann Constant

a = 7.565767e-16        #J/m^3 K^4

gff = 1
fpp = 1
f3a = 1

R_sol = 6.9634e8        #Solar Radius
P_sol = 2.65e16         #Solar Core Pressure
M_sol = 1.98847e30      #Solar Mass
L_sol = 3.916e26        #Solar Luminosity
T_sol = 15e6            #Solar Core Temperature

rho_sol = M_sol/(R_sol**3)      #Solar Density constant
eps_sol = L_sol/M_sol           
kc_sol = (R_sol*R_sol)/M_sol

#Constants for ODE's
c1 = - (G * rho_sol * M_sol )/ (P_sol * R_sol)
c2 = (4 * np.pi * R_sol * R_sol * R_sol * rho_sol) /M_sol
c3 = (4 * np.pi * R_sol * R_sol * R_sol * rho_sol * eps_sol) /L_sol
c4 = (-3 * L_sol * rho_sol * kc_sol) /(16 * np.pi * a * c * (T_sol**4) * R_sol) #Radiative Transport
c5 = - (1 - 1/gamma) * (Mh * G * M_sol) /(k * T_sol * R_sol)                    #Convective Transport