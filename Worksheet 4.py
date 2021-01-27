# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 11:34:12 2021

@author: ryank
"""

import numpy as np
import matplotlib.pyplot as plt


#Q1 sort random numbers into a set of bins
n = 100000
nbins = 10


numrand = np.random.rand(n)


plt.plot(numrand, ",", alpha = 0.1)
plt.show()

bins = np.zeros(nbins)

for i in range(len(numrand)):
    
    b = int((numrand[i]-np.min(numrand))/(np.max(numrand)+1e-10-np.min(numrand))*nbins)
    bins[b] += 1
    
    
plt.bar(range(nbins), bins)
plt.title('random')
plt.show()





#Q2 sort random gaussian distribution of numbers into a set of bins
bins = np.zeros(nbins)

nGauss = np.random.normal(0, 1, size = n)

for i in range(len(nGauss)):
    
    b = int((nGauss[i]-np.min(nGauss))/(np.max(nGauss)+1e-10-np.min(nGauss))*nbins)
    
    bins[b] += 1
    
    
plt.bar(range(nbins), bins)
plt.title('gaussian')
plt.show()





#Q3 approximate pi using a uniform distribution

points = 1000
insidePoints = 0

for i in range(0,points):
    x = np.random.rand()
    y = np.random.rand()
    
    if np.sqrt((x**2)+(y**2)) <= 1:
        insidePoints += 1
        
pi = 4*insidePoints/points

print('Our approximation of pi using', points, 'points is', pi)



















