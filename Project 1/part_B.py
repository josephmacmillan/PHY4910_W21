# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 17:00:44 2021

@author: Scott
"""
#Part B


#arrays with 20 zeros in them for the radius, temp, luminosity, mass 
r1 = np.zeros(20)
T1 = np.zeros(20)
L1 = np.zeros(20)
M1 = np.zeros(20)
p1 = np.zeros(20)

x = 0.7
x1 = 0.86
sigma1 = 5.67*(10**-5)
#loops to run randomize solutions 20 times and puts the last entries of radius, mass, pressure, Temp & Luminosity into the arrays
for a in range(int(x),int(x1)): 
    for i in range(20): #should loop 20 numbers between our X value of 0.7 to 0.86 
        r, T2, L, M, p2 = SymStatStar(X[a], Y[a], Z[a], T[i], p[i])
        t_eff = (L[-1]/(4*np.pi*r[-1]**2*sigma1))**(0.25)
        r1[i] = r[-1]
        T1[i] = t_eff
        L1[i] = L[-1]
        M1[i] = M[-1]
        p1[i] = p2[-1]

#scatter plot/HR diagram of the 20 models to see how they look 
plt.scatter(T1, L1,marker = "o", color = "k")
#plt.plot(p1, r1)
plt.xlabel('temperature')
plt.ylabel('Luminosity')
plt.show()