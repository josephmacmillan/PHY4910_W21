from nbody import *
from random import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3

from sys import argv
from time import sleep

fig = plt.figure()
ax = p3.Axes3D(fig)
#ax.set_aspect('equal')
#ax.set_xlim(-1000, 1000)
#ax.set_ylim(-1000, 1000)
ax.set_xlim3d([0,50])
ax.set_ylim3d([0,50])
ax.set_zlim3d([0,50])
ax.view_init(elev=45)

s = System.read(argv[1])
x = s.all_x()
y = s.all_y()
z = s.all_z()

line, = ax.plot(x, y, color="black", marker="o", markersize=1, linestyle='none', mfc="grey", mec="black", mew=0.5)

def update(fname):
    #sleep(0.1)
    s = System.read(fname)
    x = s.all_x()
    y = s.all_y()
    z = s.all_z()
    
    line.set_data(x, y)
    line.set_3d_properties(z)    

Writer = animation.writers['ffmpeg']
writer = Writer(fps=24, metadata=dict(artist='Me'), bitrate=-1)

dpi = 480
if len(argv) > 2:
    ani = animation.FuncAnimation(fig, update, frames=argv[2:], interval=20)


    
ani.save('willowisps.mp4', writer=writer, dpi=dpi)

plt.show()