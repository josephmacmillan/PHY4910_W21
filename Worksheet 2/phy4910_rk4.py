import numpy as np
import matplotlib.pyplot as plt

def ode_solver(x_start, x_stop, h, y0, z0, f, g):
    
    x = np.arange(x_start,x_stop,h)
    y = np.zeros(len(x))
    z = np.zeros(len(x))
    
    y[0] = y0
    z[0] = z0
    

    for i in range(0,len(x)-1):
        l1 = h*g(x[i],y[i],z[i])
        k1 = h*f(x[i],y[i],z[i])
        
        l2 = h*g(x[i] + h/2, y[i] + k1/2,z[i] + l1/2)
        k2 = h*f(x[i] + h/2, y[i] + k1/2,z[i] + l1/2)
        
        l3 = h*g(x[i] + h/2, y[i] + k2/2,z[i] + l2/2)
        k3 = h*f(x[i] + h/2, y[i] + k2/2,z[i] + l2/2)
        
        l4 = h*g(x[i] + h, y[i] + k3,z[i] + l3)
        k4 = h*f(x[i] + h, y[i] + k3,z[i] + l3)
        
        y[i+1] = y[i] + k1/6 + k2/3 + k3/3 + k4/6
        z[i+1] = z[i] + l1/6 + l2/3 + l3/3 + l4/6
    
    return x,y,z


