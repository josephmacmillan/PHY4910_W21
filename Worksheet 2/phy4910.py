# -*- coding: utf-8 -*-
import numpy as np


def odesolver(start, stop, stepsize, y0, z0, f, g):
    '''
    Ode Solver that utilizes Euler's Method and RK4.
    '''
    
    '''
    Euler's Method'
    '''
    # Initializes the x,y,z components as well as the ICs
    xn = np.arange(start, stop, stepsize)
    yn = np.zeros(len(xn))
    zn = np.zeros(len(xn))
    
    yn[0] = y0
    zn[0] = z0
    
    for i in range(len(xn) -1):
        # Loops to find new y and z based on the functions
        k = stepsize * f(xn[i], yn[i], zn[i])
        l = stepsize * g(xn[i], yn[i], zn[i])
        
        yn[i+1] = yn[i] + k
        zn[i+1] = zn[i] + l
        
    
    '''
    Runge-Kutta Method
    '''
    # Initializes the x,y,z components as well as the ICs
    xr = np.arange(start, stop, stepsize)
    yr = np.zeros(len(xr))
    zr = np.zeros(len(xr))
    
    yr[0] = y0
    zr[0] = z0
        
    for i in range(len(xr) - 1):
        # Loops to find new y and z based on functions and ki/li values
        k1 = stepsize * f(xr[i], yr[i], zr[i])
        l1 = stepsize * g(xr[i], yr[i], zr[i])
            
        k2 = stepsize * f(xr[i] + stepsize/2, yr[i] + k1/2, zr[i] + l1/2)
        l2 = stepsize * g(xr[i] + stepsize/2, yr[i] + k1/2, zr[i] + l1/2)
            
        k3 = stepsize * f(xr[i] + stepsize/2, yr[i] + k2/2, zr[i] + l2/2)
        l3 = stepsize * g(xr[i] + stepsize/2, yr[i] + k2/2, zr[i] + l2/2)
            
        k4 = stepsize * f(xr[i] + stepsize, yr[i] + k3, zr[i] + l3)
        l4 = stepsize * g(xr[i] + stepsize, yr[i] + k3, zr[i] + l3)
            
        yr[i+1] = yr[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
            
        zr[i+1] = zr[i] + (1/6)*(l1 + 2*l2 + 2*l3 + l4)
        
    return xn,yn,zn,yr,zr
