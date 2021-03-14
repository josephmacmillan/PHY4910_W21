# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 15:54:43 2021

@author: Scott
"""
import matplotlib.pyplot as plt
import numpy as np

x, y = np.loadtxt('HR_Data.txt', delimiter=',', unpack=True)
plt.scatter(x,y, label='HR diagram')
plt.gca().invert_xaxis()
plt.xlim([40000,2500])
plt.ylim([0.00001,1000000])

plt.xlabel('Temp')
plt.ylabel('Luminosity')
plt.title('HR diagram')
plt.legend()
plt.show()