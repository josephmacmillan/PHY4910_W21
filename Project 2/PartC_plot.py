from nbody import System
import numpy as np
import matplotlib.pyplot as plt


"""
Plots the positions of particles at any frame of the simulation
"""

s = System.read("./ori/Naru/Naru_1.000.dat") #Reads the data file of the final frame of the simulation for a box with 64**3 particles.
x = s.all_x()
y = s.all_y()

plt.plot(x,y, marker = 'o', linestyle='none', color='k', markersize=3)
plt.show()