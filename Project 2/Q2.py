# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 15:25:07 2021

@author: Scott
"""
from nbody3 import *
s = System.read_binary("init_z30_L50_N32.dat")

for i in range(len(s.particles)):
    for j in range(3):
        s.particles[i].position[j] = 2*(s.particles[i].position[j])
s.write("PD_Q2.dat")
#This doubled the box position by a multiple of 2      