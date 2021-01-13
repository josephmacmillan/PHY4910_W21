import numpy as np
import matplotlib.pyplot as pl

x = np.linspace(0,1,100)
y = x*np.exp(-1*x*x)

z = np.column_stack((x,y))
np.savetxt('data.txt',z)