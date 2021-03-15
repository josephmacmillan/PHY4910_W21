import matplotlib.pyplot as plt
import numpy as np
from constants import *

#Loads data from model of N stars
data = np.loadtxt("comptest.txt",unpack=True)

#Places each array into respective variables
r = data[0]
P = data[1]
M = data[2]
L = data[3]
T = data[4]


#Plots HR Diagram with N stars
fig,ax1 = plt.subplots()
ax1.set_xlabel("Temperature")
ax1.set_ylabel("Luminosity")
ax1.scatter(T*T_sol,L*L_sol)
ax1.invert_xaxis()
#ax1.set_xscale("log")
#ax1.set_yscale("log")

plt.show()

