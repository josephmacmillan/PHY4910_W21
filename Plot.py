# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 12:47:38 2021

@author: Darryen
"""
import matplotlib.pyplot as plt
import numpy as np

inputfile = input("Input the data file you want to plot: \n")
data = np.loadtxt(inputfile, unpack = True)
# We take the data from data.txt and then we make sure that the data is unpacked. 
# The returned array is transposed, so that arguments may be unpacked using x, y, z = loadtxt(...). 
# When used with a structured data-type, arrays are returned for each field.

fig1 = plt.figure()
# First we make a figure that can be saved or modified.

choiceoflimits = input("Do you want to manually modify the limits? Press y for yes, press anything else for no")

if choiceoflimits == "y":
    choicexlimitleft = input("Choose x limit(leftside): ")
    choicexlimitright = input("Choose x limit(rightside): ")

    choiceylimitleft = input("Choose y limit(leftside): ")
    choiceylimitright = input("Choose y limit(rightside): ")
    
# Above asks for limits if user wants to input them.

plt.plot(data[0], data[1])
plt.ylabel("$f(x)$")
plt.xlabel("$x$")
if choiceoflimits == "y":
    plt.xlim((choicexlimitleft, choicexlimitright))
    plt.ylim((choiceylimitleft, choiceylimitright))
plt.show()
# Plots the respective data.

choicesave = input("Press s to save, otherwise print \n")
# Asks the user whether they want to save a pdf of the plot or just print it.

if choicesave == "s":
    fig1.savefig("Plot.pdf")

# If the user types s then they will save the pdf otherwise the plot will be printed. 