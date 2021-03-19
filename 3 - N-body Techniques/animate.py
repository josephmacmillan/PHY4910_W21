from nbody import *
from random import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from sys import argv
from time import sleep

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_aspect('equal')
ax.set_xlim(-30, 30)
ax.set_ylim(-30, 30)

s = System.read(argv[1])
x = s.all_x()
y = s.all_y()

line, = ax.plot(x, y, color="black", marker="o", markersize=5, linestyle='none', mfc="grey", mec="black", mew=0.5)

def update(fname):
    sleep(0.01)
    s = System.read(fname)
    x = s.all_x()
    y = s.all_y()
    
    line.set_data(x, y)    

if len(argv) > 2:
    ani = animation.FuncAnimation(fig, update, frames=argv[2:], interval=20)

plt.show()

