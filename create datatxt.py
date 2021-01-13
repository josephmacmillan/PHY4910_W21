

"""
PHY4910
Team Planck
Worksheet #1
"""
import numpy as np          #import library to use

x = np.linspace(0,1.0,100)  #create an array of 100 vlues from 0 to 1
y = x*np.exp(-x**2)         #create array which is a function of x

data = np.column_stack((x,y)) #stack the arrays as columns
np.savetxt('data.txt',data) #save the two-columned array as a txt file