# -*- coding: utf-8 -*-
"""
Team Plank
Worksheet 3
"""


import phy4910
import numpy as np
import matplotlib.pyplot as plt

#define value of n as given
n = 1.5

#define functions in our problem in terms of eta,rho, and sigma
def f(eta,rho,sigma):
    
    return sigma

def g(eta,rho,sigma):
    
    return ((-2*sigma)/eta) - rho**n

#call solveODE to solve for values of eta, rho, sigma
eta,rho,sigma = phy4910.solveODE(0.001,10,0.001,1,0,f,g)

#trim the part of the arrays that's not a number
condition = rho > 0.0
eta = eta[condition]
rho = rho[condition]
sigma = sigma[condition]

#print the value of the edge of the star
print("Edge of the star is at eta = ", eta[-1])
#plot the variables
plt.plot(eta,rho, color = "blue")
plt.plot(eta, sigma, color = "orange")
plt.show()

#calculate dimensionless mass
y = rho**n * eta**2
m = np.trapz(y,eta)
print("Dimensionless mass m = ", m)

# define global variables for out constants
k = 3.166 * 10**(12)
G = 6.6743 * 10**(-8)
rho_c = 4.045 * 10**(6)

lamb = (((n + 1) * k * rho_c**((1-n)/n))/ (4 * np.pi * G))**(0.5) #in cm
lamb = lamb * 10**(-5)  #convert to km
print("lambda = ", lamb)


#find real density and radius of the star

eta = eta * lamb  # km
rho = rho * rho_c
print("Radius of the star = ", eta[-1], " Km")
#plot the density vs radius with units
plt.plot(eta,rho, color = "pink")
plt.show()


#Calculate mass M
#within the formula multiply lambda by 10^5 to convert to cm
# our final asnwer for M will be in grams, so we will multiply the entire 
#formula by 10^-3 to convert to kg

M = 4 * np.pi * rho_c * (lamb * 10**(5))**3 * m * 10**(-3)
print("Mass of the star = ", M, " Kg")




