from math import cos, sin, pi, sqrt
from random import random
from nbody import *

# units in earth masses, AU, years

N = 9
t = 0.0
G = 1.184e-4
masses = [332979, 0.05528, 0.816, 1.0, 0.108, 317.9, 95.14, 14.54, 17.08]
radii = [0.0, 0.387, 0.7233, 1.0, 1.523, 5.20, 9.58, 19.20, 30.05]
velocity = [0.0, 9.994, 7.38, 6.283, 5.08, 2.76, 2.046, 1.43, 1.14]

Px = 0.0
Py = 0.0
Pz = 0.0
mtot = 0.0
parts = []
for i in range(N):
    angle = 2.0 * pi *random()
    x = radii[i] * cos(angle)
    y = radii[i] * sin(angle)
    z = 0
    
    vx = -velocity[i] * sin(angle)
    vy = velocity[i] * cos(angle)
    vz = 0.0
    
    m = masses[i]
    Px += m * vx
    Py += m * vy
    Pz += m * vz
    mtot += m
    
    parts.append(Particle(m, [x, y, z], [vx, vy, vz]))
    
s = System(parts, t, G)

# what is the total momentum of the system?
P = sqrt(Px**2 + Py**2 + Pz**2)
print(P)

# that's no good, switch to centre of momentum frame
Vx = Px / mtot
Vy = Py / mtot
Vz = Pz / mtot

for p in s.particles:
    p.velocity[0] -= Vx
    p.velocity[1] -= Vy
    p.velocity[2] -= Vz
    
s.write("solar_system.dat")

Px = 0
Py = 0
Pz = 0
xcm = 0
ycm = 0
zcm = 0
for p in s.particles:
    Px += p.mass * p.vx()
    Py += p.mass * p.vy()
    Pz += p.mass * p.vz()

    xcm += p.mass * p.x()
    ycm += p.mass * p.y()
    zcm += p.mass * p.z()

P = sqrt(Px**2 + Py**2 + Pz**2)
xcm /= mtot
ycm /= mtot
zcm /= mtot
Rcm = sqrt(xcm**2 + ycm**2 + zcm**2)
print(P)
print(xcm, ycm, zcm, Rcm)
