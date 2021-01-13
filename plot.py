# -*- coding: utf-8 -*-
"""
PHY4910
Team Planck
Worksheet #1
"""

#import libraries needed
import numpy as np
import matplotlib.pyplot as plt


#load saved txt file
data = np.loadtxt('data.txt', unpack=True)

#extract each column as an array
x= data[0]
y= data[1]

#plot the arrays
plt.plot(x,y)
#set axis limits
plt.xlim(-0.1,1.2*x[-1])
plt.ylim(-0.1,0.8)

#label axes and make a title for the plot
plt.xlabel('x')
plt.ylabel('f(x)',fontsize=16)
plt.title('f(x) vs x')
#display the plot
plt.show()
