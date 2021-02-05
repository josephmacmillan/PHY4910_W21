# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 11:21:57 2021

@author: ryank
"""

import numpy as np
import matplotlib.pyplot as plt

#worksheet 5, part A

rng = np.random.default_rng()

#returns random theta and phi values 
def pick_direction():
    
    #correction for the distribution to prevent clustering at the peaks 
    #see link: http://corysimon.github.io/articles/uniformdistn-on-sphere/
    theta = np.arccos(1 - 2 * rng.random())

    theta = np.pi*rng.random() 
    phi = 2*np.pi*rng.random() 
    
    return theta, phi

print(pick_direction())


#Tests pick_direction()
def test_dir(n):
    arrtheta = np.zeros(n)
    arrphi = np.zeros(n)
    
    for i in range (len(arrtheta)):
        arrtheta[i], arrphi[i] = pick_direction()
        
        
    #Convert to cartesian coordinates, setting r = 1
    x = np.sin(arrtheta)*np.cos(arrphi)
    y = np.sin(arrtheta)*np.sin(arrphi)
    z = np.cos(arrtheta)    
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x,y,z,marker ='.')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()
    
test_dir(10000)
    
    
#worksheet 5, Part B

#returns random optical depth distance (tau)
def pick_optical_depth():
    x = rng.random()
    tau = -np.log(1-x)
    
    return tau


#Tests pick_optical_depth() 
def test_op_depth(n,nbin):
    arrdepth = np.zeros(n)
    
    #store random optical depth distance values in arrdepth 
    for i in range(n):
        arrdepth[i] = pick_optical_depth()
    
    dmax = np.amax(arrdepth)
    dmin = np.amin(arrdepth)
    
    bins = np.zeros(nbin)
    
    #bin the optical depth data
    for i in range(n):
        b = int((arrdepth[i] - dmin)/(dmax + 1e-7 - dmin)*nbin)
        bins[b] += 1
        
    plt.bar(range(nbin), bins)
    plt.xlabel('bins')
    plt.ylabel('optical depth distance')
    plt.show()
    
test_op_depth(10000,20)
    
    
    
    
#worksheet 5, Part C


def move_photon(tmax, zmax):
    x = []
    y = []
    z = []
    
    #starting at the origin
    x.append(0)
    y.append(0)
    z.append(0)
    
    #loop until photon moves from origin to the surface
    while True: 
        theta, phi = pick_direction()
        opdepth = pick_optical_depth()
        
        
        #how far the photon moves in the atmophere
        s = opdepth/tmax
        
        #how far the photon has moved since last time (not current position in atmosphere)
        dx = s * np.sin(theta)*np.cos(phi)
        dy = s * np.sin(theta)*np.sin(phi)
        dz = s * np.cos(theta)
    
        #returns last scattering direction of the photon at the surface
        x.append(x[-1] + dx)
        y.append(y[-1] + dy)
        z.append(z[-1] + dz)
        
        #fancy fstring stuff
        print(f"Photon Position{x[-1]},{y[-1]},{z[-1]}")
        
        #if photon happens to go in the wrong direction, reset position of the photon
        if z[-1] < 0:
            x.clear()
            y.clear()
            z.clear()
            x.append(0)
            y.append(0)
            z.append(0)
        
        #once the photon hits the surface, return the data
        if z[-1] > zmax:
            return x,y,z,theta,phi
    
    
x,y,z,theta,phi = move_photon(10,1)    

#plot path of the photon from origin to the surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x,y,z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
    
    
    
    
    
    