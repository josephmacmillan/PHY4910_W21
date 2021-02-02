# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 11:34:12 2021

"""

import numpy as np
import matplotlib.pyplot as plt


#Q1 sort random numbers into a set of bins
n = 100000
nbins = 10


#create a uniform distribution of random numbers between 0 and 1
numrand = np.random.rand(n)


#plot the distribution of random numbers (alpha lowers transparency of datapoints)
plt.plot(numrand, ",", alpha = 0.1)
plt.show()


#create an array filled with zeros with length nbins
bins = np.zeros(nbins)


#Sorts values from the random uniform distribution into bins
#loop calculates b for each value of numrand and adds 1 to the bin it falls into (see lecture 3 for b equation)
for i in range(len(numrand)):
    
    b = int((numrand[i]-np.min(numrand))/(np.max(numrand)+1e-10-np.min(numrand))*nbins)
    bins[b] += 1
    

#plot bar graph to compare distribution among bins    
plt.bar(range(nbins), bins)
plt.title('random')
plt.ylabel('value of random number')
plt.xlabel('bins')
plt.show()





#Q2 sort random gaussian distribution of numbers into a set of bins

#create an array filled with zeros with length nbins
bins = np.zeros(nbins)


#create a uniform gaussian (normal) distribution of length n
nGauss = np.random.normal(0, 1, size = n)


#Sorts values from the random uniform distribution into bins
#loop calculates b for each value of nGauss and adds 1 to the bin it falls into
for i in range(len(nGauss)):
    
    b = int((nGauss[i]-np.min(nGauss))/(np.max(nGauss)+1e-10-np.min(nGauss))*nbins)
    bins[b] += 1
    
    
#plot bar graph to compare distribution among bins     
plt.bar(range(nbins), bins)
plt.title('gaussian')
plt.show()





#Q3 approximate pi using a uniform distribution of random numbers
points = 1000
insidePoints = 0

#imagine a circle of radius r = 1 enclosed by a square with sidelengths 2r

#loop randomly distributes points (values between 0 and 1) and calculates pi using... 
#...the numbers of points that land within the unit circle (Monte Carlo Method)
for i in range(0, points):
    
    #generates random coordinate points with x and y values between 0 and 1
    x = np.random.rand()
    y = np.random.rand()
    
    #if the line drawn out between the x and y coordinate and the origin is less than the unit radius...
    #...the point lands within the circle. The number of inside points is incremented by 1 
    if np.sqrt((x**2)+(y**2)) <= 1:
        insidePoints += 1
        

#dividing the expression for the area of a circle by the expression for the area of a circle gives Pi/4
#Same ratio can be used between the # of points that land in the circle and the # of points that land in the square
pi = 4*(insidePoints/points)

print('Our approximation of pi using', points, 'points is', pi)



















