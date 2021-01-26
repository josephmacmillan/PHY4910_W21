# -*- coding: utf-8 -*-
"""
Tech of modern astro
team Planck
Worksheet 2

Author: Georges Karagozian
"""

import numpy as np


def solveODE(x_start, x_stop, h, y0, z0, f, g, Type = "rk"):
   
    
    #initialize the arrays for each variable
    x = np.arange(x_start,x_stop,h)
    N = len(x)
    y = np.zeros(N)
    z = np.zeros(N)
    
    #set initial conditions
    y[0] = y0
    z[0] = z0
    
    # if conditional statement to decide which solution method to use
    
    # default to runge-kutta
    if Type == "rk" :
        for i in range(0,N-1):
            # Loops to find new y and z based on functions and ki/li values
            k1 = h * f(x[i],y[i],z[i])
            l1 = h * g(x[i],y[i],z[i])
            
            k2 = h * f(x[i] + h/2, y[i] + k1/2, z[i] + l1/2)
            l2 = h * g(x[i] + h/2, y[i] + k1/2, z[i] + l1/2)
            
            k3 = h * f(x[i] + h/2, y[i] + k2/2, z[i] + l2/2)
            l3 = h * g(x[i] + h/2, y[i] + k2/2, z[i] + l2/2)
            
            k4 = h * f(x[i] +h, y[i] +k3, z[i] +l3)
            l4 = h * g(x[i] +h, y[i] +k3, z[i] +l3)
            
            
            #use kn's and ln's (derivatives) to compute yn+1 and zn+1 (following points)
            y[i+1] = y[i] + k1/6 + k2/3 + k3/3 + k4/6
            z[i+1] = z[i] + l1/6 + l2/3 + l3/3 + l4/6
    
    # euler method
    elif Type == "eu" :
        
        for i in range(0, N-1):
            
            y[i+1] = y[i] + z[i] * h
            z[i+1] = z[i] - y[i] * h
    
    # in case a random string is passed display error    
    elif Type != "eu" and Type != "rk":
        print("Error! Please enter a valid type of solution method.")
    
    # return arrays of x,y,z variables
    return x,y,z    
        
        
        
def plot(x,y, xlabel, ylabel, save = False, name = " plot.pdf"):
    import matplotlib.pyplot as plt
    
    #plot the arrays
    plt.plot(x,y)
    

    #label axes and make a title for the plot
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    #save plot if indicated in arguments
    if save == True:
        plt.savefig(name)
    #display the plot
    plt.show()
        
    return None    
        
 


def Nonrel_WhiteDwarf(x_start, x_stop, h, y0, z0, f, g, n, k, G, rho_c, lamb):
    #import libraries
    import matplotlibe.pyplot as plt
    import numpy as np
    
    #call solveODE to solve for values of eta, rho, sigma
    eta,rho,sigma = solveODE(x_start, x_stop, h, y0, z0,f,g)

    #trim the part of the arrays that's not a number
    condition = rho > 0.0
    eta = eta[condition]
    rho = rho[condition]
    sigma = sigma[condition]
    
    #print the value of the edge of the star
    print("Edge of the star is at eta = ", eta[-1])
    #plot the variables
    plt.plot(eta,rho, color = "blue", label = '$\\rho$')
    plt.plot(eta, sigma, color = "orange", label = "$\sigma$")
    plt.legend()
    plt.xlabel('$\eta$')
    plt.show()
    
    #calculate dimensionless mass
    y = rho**n * eta**2
    m = np.trapz(y,eta)
    print("Dimensionless mass m = ", m)
    

    #write formula for lambda then convert to km
    lamb = (((n + 1) * k * rho_c**((1-n)/n))/ (4 * np.pi * G))**(0.5) #in cm
    lamb = lamb * 10**(-5)  #convert to km
    print("lambda = ", lamb)
    
    
    #find real density and radius of the star
    
    eta = eta * lamb  # convert to km using the constant lambda
    rho = rho * rho_c 
    print("Radius of the star = ", eta[-1], " Km")
    
    #plot the density vs radius with units
    plt.plot(eta,rho, color = "pink", label = "Real eta vs. rho")
    plt.legend()
    plt.show()
    
    
    #Calculate mass M
    #within the formula multiply lambda by 10^5 to convert to cm
    # our final asnwer for M will be in grams, so we will multiply the entire 
    #formula by 10^-3 to convert to kg
    
    M = 4 * np.pi * rho_c * (lamb * 10**(5))**3 * m * 10**(-3)
    print("Mass of the star = ", M, " Kg")
        
        
        
    return None
