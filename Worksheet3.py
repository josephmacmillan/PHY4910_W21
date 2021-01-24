# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 11:22:59 2021

@author: Team Chandra
"""
import numpy as np
import matplotlib.pyplot as plt 
import phy4910 as ph

def f(x,y,z):
    return z

def g(x,y,z):
    n = 1.5
    return (-2/x) * z - y**n 

# Use odesolver to find surface of the star

x,y, yr = ph.odesolver(0.000001, 4, 0.001, 1, 0, f,g)

condition = yr > 0.0
x = x[condition]
y = y[condition]
yr = yr[condition]


print("\nThe surface of the star is", x[-1])

plt.plot(x,yr)
plt.xlabel('$\eta$')
plt.ylabel('$\\rho$')
plt.show()

# Part 2)
# Use trapz to numerically integrate to find dimensionless mass
n = 1.5
m = np.trapz(yr**n * x**2, x)
print("\nThe dimensionless mass is", m)


#Part 3)
# Initialize the constants
k = 3.166E12
pc = 4.045E6
G = 6.6743E-8
lam = (((n+1) * k * np.power(pc, (1 - n) /n))/(4*np.pi*G))**(1/2)
lam *= 10**(-5)
print("\nThis is the radial scale factor", lam, "km")

x *= lam
yr *= pc
plt.plot(x,yr)
plt.xlabel('$r$')
plt.ylabel('$\\rho$')
plt.show()

print("\nThe radius of the white dwarf is", x[-1], "km")

# Change rho to fit the units to get kg value
M = 4*np.pi*pc*(lam * 10**(5))**3 * m * 10**(-3)

print("\nThe mass of the White Dwarf is", M, "kg")
print("\nThis is a decent approximation because the Chandrasekhar limit is", 2.7E30, "kg")

print("\n \n Part B)")

# B
#--------------------------------------------------------------------------------

# Part B is much of the same thing
def f1(x,y,z):
    return z

def g1(x,y,z):
    n = 3
    return (-2/x) * z - y**n 


x1,y1, yr1 = ph.odesolver(0.000001, 4, 0.001, 1, 0, f1,g1)

condition = yr1 > 0.0
x1 = x1[condition]
y1 = y1[condition]
yr1 = yr1[condition]


print("\nThe surface is", x1[-1])

plt.plot(x1,yr1)
plt.xlabel('$\eta$')
plt.ylabel('$\\rho$')
plt.show()

# Part 2)
n1 = 3
m1 = np.trapz(yr1**n1 * x1**2, x1)
print("\nThe dimensionless mass is", m1)


#Part 3)
k1 = 4.936E14
pc1 = 53.31E6
lam1 = (((n1+1) * k1 * np.power(pc1, (1 - n1) /n1))/(4*np.pi*G))**(1/2)
lam1 *= 10**(-5)
print("\nThis is the radial scale factor", lam1, "km")

x1 *= lam1
yr1 *= pc1
plt.plot(x1,yr1)
plt.xlabel('$r$')
plt.ylabel('$\\rho$')
plt.show()

print("\nThe radius of the white dwarf is", x1[-1], "km")


M1 = 4*np.pi*pc1*(lam1*10**(5))**3 * m1 *10**(-3)

print("\nThe mass of the White Dwarf is", M1, "kg")
print("\nThis is a worst approximation because of the Chandrasehkar limit.")
