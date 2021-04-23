import random as rand
from nbody import Particle
from nbody import System
import numpy as np

"""
Produces two Initial Conditions data files for particles with uniform and random positions
"""

L = 1                   #Box length
Ng = 64                 #Number of grid points

Np = 32                 #Number of particles in the box
M = 1                   #Total mass of particles
m = M/(Np**3)           #Mass of individual particles

#Uniform particle disturbution
order = []              #List to append Particle attributes

for i in range(Ng):     #assigns uniform particle positions to particles
    for j in range(Ng):
        for k in range(Ng):
            x = i
            y = j
            z = k
            particle = Particle(m, [x,y,z], [0,0,0])                         
            order.append(particle)
System(order).write('order.dat')    #writes the System data of uniform particle system 
            


#Random particle distribution
chaos = []              #List to append Particle Attributes
for i in range(Ng):     #assigns random particle positions to particles
    for j in range(Ng):
        for k in range(Ng):
            x = rand.random()*Ng
            y = rand.random()*Ng
            z = rand.random()*Ng
            particle = Particle(m, [x,y,z], [0,0,0])                         
            chaos.append(particle)
            
ch = System(chaos)
ch.write('chaos.dat')           #writes the System data of random particle system
