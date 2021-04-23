import random as rand
from nbody import Particle
from nbody import System
import numpy as np

"""
Reads data from initial conditions file and scale positions twice as normal for a box with doubled length.
Write new particle IC's to new file'
"""

s = System.read_binary('init_z30_L50_N32.dat')

for i in range(len(s.particles)):
    for j in range(3):
        s.particles[i].position[j] += s.particles[i].position[j]
    
s.write('new.dat')
