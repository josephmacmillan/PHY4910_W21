# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 01:49:42 2021

@author: Darryen
"""
import Mesh
import numpy as np
import matplotlib.pyplot as plt
import nbody 

# Case 1 

# Need to make the system of evenly distributed particles on a 3d grid

particles = []


# 16^3 particles 
N = 16

# Side length of one
L = 2

# A loop that evenly distributes particles in a 3D grid
for i in range(N):
    x = i / (N*L)
    for j in range(N):
        y = j / (N*L)
        for k in range(N):
            z = k /(N*L)
            # Places a particle with masses that total to 1, and places them with correct positions
            particles.append(nbody.Particle(1/(N**3), [x,y,z], [0,0,0]))
            


# Creates the mesh that we will evolve our system in
mesh = Mesh.Mesh(L, N*2)

# Creates our system of particles 
system = nbody.System(particles, t = 0, G = 1)

a = 1/31.0 # z = 30 # Current age of the universe
da = 0.01 # Acts as our universe timestep

# Evolves our system by calculating the accelerations and updating the particles
while a <= 1.0:
    accels = mesh.calc_accels(system, a)
    mesh.move_particles(system, accels, a, da)
    
    
    system.write(f"output_{a:5.3f}")
    a += da
    
    '''
# Case 2


# This time we have a random ensemble of particles

particles = []

N = 16

# A loop that randomly distributes particles in a 3D grid
for i in range(N):
    x = np.random.random() * L
    y = np.random.random() * L
    z = np.random.random() * L
    particles.append(nbody.Particle(1/(N**3), [x,y,z], [0,0,0]))
            


# Creates the mesh that we will evolve our system in
mesh = Mesh.Mesh(L, N*2)

# Creates our system of particles 
system = nbody.System(particles, t = 0, G = 1)

a = 1/31.0 # z = 30 # Current age of the universe
da = 0.01 # Acts as our universe timestep

# Evolves our system by calculating the accelerations and updating the particles
while a <= 1.0:
    accels = mesh.calc_accels(system, a)
    mesh.move_particles(system, accels, a, da)
    
    
    system.write(f"output_{a:2.3f}")
    a += da
    
'''
