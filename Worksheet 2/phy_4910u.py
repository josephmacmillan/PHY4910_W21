# -*- coding: utf-8 -*-
"""
Spyder Editor

Author: Ryan Kwok
"""

import matplotlib.pyplot as plt
import numpy as np

def ode_solver_euler(x_start, x_stop, h):


    #make an array of x-values with certain step size 
    #h = 0.01
    x = np.arange(x_start, x_stop, h)
    
    
    #make arrays of zeros to be filled for y and z
    y = np.zeros(len(x))
    z = np.zeros(len(x))
    
    
    #define functions for each of the first order ODEs
    def f(x,y,z):
        return z
    
    def g(x,y,z):
        return -y
    
    
    #boundary conditions at xn = 0
    y[0] = 1
    z[0] = 0
    
    
    #use loop to fill in the array of zeros for y and z
    for i in range(len(x)-1):
        kn = h*f(x[i],y[i],z[i])
        ln = h*g(x[i],y[i],z[i])
        
        y[i+1] = y[i] + kn
        z[i+1] = z[i] + ln
    
    #return values
    return x, y, z



x,y,z = ode_solver_euler(0, 50, 0.01)
plt.title('euler')
plt.plot(x,y)
plt.show()




def ode_solver_rungeKutta(x_start, x_stop, h):
    #make an array of x-values with certain step size 
    #h = 0.01
    x = np.arange(x_start, x_stop, h)
    
    
    #make arrays of zeros to be filled for y and z
    y = np.zeros(len(x))
    z = np.zeros(len(x))
    
    
    #define functions for each of the first order ODEs
    def f(x,y,z):
        return z
    
    def g(x,y,z):
        return -y
    
    #boundary conditions at xn = 0
    y[0] = 1
    z[0] = 0
    
    #use loop to compute values for y and z using the runge kutta method
    #one derivative at the start, two at midpoint, one at the end
    
    for i in range(len(x)-1):
        #compute all kn's and ln's (derivatives)
        
        #one derivative at the start
        k1 = h*f(x[i],y[i],z[i])
        l1 = h*g(x[i],y[i],z[i])
        
        #two derivatives at the midpoint
        k2 = h*f(x[i]+(h/2), y[i]+(k1/2), z[i]+(l1/2))
        l2 = h*g(x[i]+(h/2), y[i]+(k1/2), z[i]+(l1/2))
        
        k3 = h*f(x[i]+(h/2), y[i]+(k2/2), z[i]+(l2/2))
        l3 = h*g(x[i]+(h/2), y[i]+(k2/2), z[i]+(l2/2))
        
        #one derivative at the endpoint
        k4 = h*f(x[i]+h, y[i]+k3, z[i]+l3)
        l4 = h*g(x[i]+h, y[i]+k3, z[i]+l3)
        
        
        
        #use kn's and ln's (derivatives) to compute yn+1 and zn+1 (following points)
        y[i+1] = y[i] + (k1/6) + (k2/3) + (k3/3) + (k4/6)
        z[i+1] = z[i] + (l1/6) + (l2/3) + (l3/3) + (l4/6)
        
        
    #return values
    return x, y, z


x,y,z = ode_solver_rungeKutta(0, 50, 0.01)
plt.title('rungekutta')
plt.plot(x,y)
plt.show()

























