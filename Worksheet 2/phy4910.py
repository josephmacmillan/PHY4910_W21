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
        
        
        
        
        
        
        
