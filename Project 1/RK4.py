# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 20:58:46 2021

@author: Moriah
"""
import numpy as np 
 

""" 

Solves a coupled system of 4 ODEs using runge kutta.
 
Takes as arguments:
  x_start - starting point for independent coordinate
  x_end - ending point
  h - difference between x_i and x_i+1 (i.e., delat x)
  y1_0 - initial value for first variable y1(x_start)
  y2_0 - initial value for second variable y2(x_start)
  y3_0 - initial value for third variable y3(x_start) 
  y4_0 - initial value for fourth variable y4(x_start)
  b - function for derivative of first variable (i.e., b = dy1/dx)
  c - function for derivative of second vaiable (i.e., c = dy2/dx)
  f - function for derivative of third variable (i.e., f = dy3/dx)
  g - function for deriviative of fourth variable (i.e., g = dy4/dx)
  
  
  returns five arrays, x[0,N-1], y1[0,N-1], y2[0,N-1], y3[0,N-1], and y4[0,N-1].
"""


def ode_rk4(x_start, x_end, h, y1_0, y2_0, y3_0, y4_0, b, c, f, g, epsi):
    
    #initialize arrays 
    x = np.arange(x_start, x_end, h)
    N = len(x)
    
    y1 = np.zeros(N)
    y2 = np.zeros(N)
    y3 = np.zeros(N)
    y4 = np.zeros(N)
    
    
    #Boudary Conidtions
    y1[0] = y1_0 #Preasure
    y2[0] = y2_0 #Mass
    y3[0] = y3_0 #Luminosity
    y4[0] = y4_0 #Temperature
    
    for i in range(0, N-1):
        k1 = h * b(x[i], y1[i], y2[i], y3[i], y4[i])
        l1 = h * c(x[i], y1[i], y2[i], y3[i], y4[i])
        j1 = h * f(x[i], y1[i], y2[i], y3[i], y4[i])
        q1 = h * g(x[i], y1[i], y2[i], y3[i], y4[i])
        
        
        k2 = h * b(x[i] + 0.5 * h, y1[i] + 0.5 * k1, y2[i] + 0.5 * l1, y3[i] + 0.5 *j1, y4[i] + 0.5*q1)
        l2 = h * c(x[i] + 0.5 * h, y1[i] + 0.5 * k1, y2[i] + 0.5 * l1, y3[i] + 0.5 *j1, y4[i] + 0.5*q1)
        j2 = h * f(x[i] + 0.5 * h, y1[i] + 0.5 * k1, y2[i] + 0.5 * l1, y3[i] + 0.5 *j1, y4[i] + 0.5*q1)
        q2 = h * g(x[i] + 0.5 * h, y1[i] + 0.5 * k1, y2[i] + 0.5 * l1, y3[i] + 0.5 *j1, y4[i] + 0.5*q1)
        
        
        k3 = h * b(x[i] + 0.5 * h, y1[i] + 0.5 * k2, y2[i] + 0.5 * l2, y3[i] + 0.5 * j2, y4[i] + 0.5 *q2)
        l3 = h * c(x[i] + 0.5 * h, y1[i] + 0.5 * k2, y2[i] + 0.5 * l2, y3[i] + 0.5 * j2, y4[i] + 0.5 *q2)
        j3 = h * f(x[i] + 0.5 * h, y1[i] + 0.5 * k2, y2[i] + 0.5 * l2, y3[i] + 0.5 * j2, y4[i] + 0.5 *q2)
        q3 = h * g(x[i] + 0.5 * h, y1[i] + 0.5 * k2, y2[i] + 0.5 * l2, y3[i] + 0.5 * j2, y4[i] + 0.5 *q2)
        
        k4 = h * b(x[i] + h, y1[i] + k3, y2[i] + l3, y3[i] + j3, y4[i] + q3)
        l4 = h * c(x[i] + h, y1[i] + k3, y2[i] + l3, y3[i] + j3, y4[i] + q3)
        j4 = h * f(x[i] + h, y1[i] + k3, y2[i] + l3, y3[i] + j3, y4[i] + q3)
        q4 = h * g(x[i] + h, y1[i] + k3, y2[i] + l3, y3[i] + j3, y4[i] + q3)
        
    
        y1[i+1] = y1[i] + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        y2[i+1] = y2[i] + (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0
        y3[i+1] = y3[i] + (j1 + 2.0 * j2 + 2.0 * j3 + j4) / 6.0
        y4[i+1] = y4[i] + (q1 + 2.0 * q2 + 2.0 * q3 + q4) / 6.0
    
        
        if y1[i+1] <= epsi and y1[i+1] <= epsi:
            print("Break", y1[i+1])
            break
        
            
     
    #returns arrays for x, y1, etc...    
    return x, y1, y2, y3, y4