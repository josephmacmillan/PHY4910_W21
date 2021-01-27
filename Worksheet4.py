# -*- coding: utf-8 -*-
"""
Astro
Worksheet 4
"""

import numpy as np
import matplotlib.pyplot as plt

N = 1000000
arr = np.random.rand(N)
arrmax = np.amax(arr) + 1e-10
arrmin = np.amin(arr)


#plt.plot(arr, range(N),",", alpha = 0.4)

N_bin = 1000

bins = np.zeros(N_bin)

for i in range(N):
    b = int(((arr[i] - arrmin)/(arrmax - arrmin)) * N_bin)
    
    bins[b] += 1
    

plt.bar(range(N_bin), bins)
plt.show()

# Part B

arr2 = np.random.normal(loc = 0, scale = 1, size = N)

arr2max = np.amax(arr2) + 1e-10
arr2min = np.amin(arr2)


plt.plot(arr2, range(N),",", alpha = 0.4)
plt.show()
N2_bin = 1000

bins2 = np.zeros(N2_bin)

for i in range(N):
    b = int(((arr2[i] - arr2min)/(arr2max - arr2min)) * N2_bin)
    
    bins2[b] += 1
    

plt.bar(range(N2_bin), bins2)
plt.show()

#part C
x = np.random.rand(N)
y = np.random.rand(N)


r = x**2 + y**2
a = np.sum(r<1.0)
pi = a/ N  * 4
print(pi)
