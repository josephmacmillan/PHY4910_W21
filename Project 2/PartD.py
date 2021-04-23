from nbody import System
import numpy as np
from PartA import *
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

"""
Plots, hightlight and prints the mass of the largest cluster found for the system with initial conditions for init_z30_L50_N32.dat
"""

fig = plt.figure()
ax = p3.Axes3D(fig)
#ax.set_aspect('equal')
#ax.set_xlim(-1000, 1000)
#ax.set_ylim(-1000, 1000)
ax.set_xlim3d([0,50])
ax.set_ylim3d([0,50])
ax.set_zlim3d([0,50])
ax.view_init(elev=45)
ax.set_ylabel('y')
ax.set_xlabel('x')
ax.set_zlabel('z')

rx = [24,33]
ry = [19,24]
rz = [27,36]
#Highlight manually selected region of cluster
for s, e in combinations(np.array(list(product(rx, ry, rz))), 2):
    if np.sum(np.abs(s-e)) == ry[1]-ry[0] or np.sum(np.abs(s-e)) == rz[1]-rz[0]:
        ax.plot3D(*zip(s, e), lw=10, color="fuchsia")

#Reads positions for data file
s = System.read('./ori/wisps/wisps_0.999.dat')
x = s.all_x()
y = s.all_y()
z = s.all_z()

#Plots position of particles in 3D
line, = ax.plot(x, y, color="black", marker="o", markersize=2, linestyle='none', mfc="k", mec="k", mew=0.5)
line.set_3d_properties(z)

plt.show()



#Counts the number of particle with in the region identified
count = 0
for i in range(len(x)):
    if 24 <= x[i] <= 33 and 19 <= y[i] <= 24 and 27 <= z[i] <= 36:
        count += 1
        
print(f'Number of particles in largest Structure = {count}')
print(f"Fraction of all particles in structure = {count/len(x)}")

#Get the Mass of one particle
mass = s.particles[0].mass

print(f"Mass of a single particle = {mass} M_sol")

#Calculate the total mass of the highlighted cluster
M_cl = count * mass

print(f"Total Mass of Cluster = {M_cl} M_sol")
