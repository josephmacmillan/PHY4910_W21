#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import argparse
from sys import argv

rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':14})
plt.rc('axes', labelsize=16)
plt.rcParams.update({'figure.autolayout': True})

parser = argparse.ArgumentParser(description='Plot some data.')
parser.add_argument('filename', type=str, help='file name containing data in columns')
parser.add_argument('-f', dest='out_file', default=None, help='Save as a file')
parser.add_argument('-x', dest='x_col', default=1, type=int, help='Column to use as x-axis')
parser.add_argument('-y', dest='y_col', default=2, type=int, help='Column to use as y-axis')
parser.add_argument('-s', dest='style', default='-', type=str, help='MatPlotLib format string for the plot)')
parser.add_argument('-l', dest='limits', type=float, default=None, nargs=4, help='Limits for the plot')
parser.add_argument('-a', dest='plot_all', action='store_true', help='Plot all columns at once')
args = parser.parse_args()


data = np.loadtxt(args.filename, unpack=True)


fig = plt.figure()
ax = fig.add_subplot(1,1,1)

if args.limits:
    ax.set_xlim(args.limits[0], args.limits[1])
    ax.set_ylim(args.limits[2], args.limits[3])
   
if not args.plot_all:
    ax.plot(data[args.x_col-1], data[args.y_col-1], args.style, color="black")
else:
    for i in range(1, len(data)):
        ax.plot(data[0], data[i], args.style)

if args.out_file:
    plt.savefig(args.out_file)    
else:
    plt.show()

