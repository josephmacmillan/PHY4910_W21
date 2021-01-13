# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 12:43:01 2021

@author: Darryen
"""
import numpy as np

x = np.linspace(0,1.0,100)

def f(x):
    y = x*np.exp(-x**2)
    return y

np.savetxt("data.txt", np.column_stack((x,f(x))))