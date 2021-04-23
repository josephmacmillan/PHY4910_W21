from nbody import *
from random import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from sys import argv
from time import sleep

"""
Reads data files of position of particles and creates a 2D animation of this.
"""

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_aspect('equal')
ax.set_xlim(0, 50)
ax.set_ylim(0, 50)

s = System.read(argv[1])        #Read data files
x = s.all_x()                   #List of x positions of all points
y = s.all_y()                   #List of y position of all points

line, = ax.plot(x, y, color="black", marker="o", markersize=1, linestyle='none', mfc="grey", mec="black", mew=0.5)

def update(fname):
    sleep(0.5)
    s = System.read(fname)
    x = s.all_x()
    y = s.all_y()
    line.set_data(x, y)    

if len(argv) > 2:
    ani = animation.FuncAnimation(fig, update, frames=argv[2:], interval=20)
    
Writer = animation.writers['ffmpeg']
writer = Writer(fps=24, metadata=dict(artist='Me'), bitrate=-1)

dpi = 480
ani.save('wisps.mp4', writer=writer,dpi= dpi)

plt.show()

