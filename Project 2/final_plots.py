# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 13:43:36 2021

@author: Scott
"""
from nbody3 import *
import matplotlib.pyplot as plt
s = System.read("Output_0.712.dat")
x = s.all_x()
y = s.all_y()
plt.plot(x, y, linestyle = "none", marker = "o")
plt.show()