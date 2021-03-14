import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pandas_ods_reader import read_ods

path = "/home/ryank/Documents/mesadata/MESA_DATA_TABLE"


# load a sheet based on its name
sheet_name = "MESA_DATA_TABLE"
df = read_ods(path, sheet_name)

# a scatter plot comparing effective Temperature and Luminosity (both log)
df.plot(kind='scatter',x='log_Teff' ,y='log_L',color='blue')

plt.title('MESA Simulation of 20 Stars with Various Masses')
plt.xlim(4.75,3.4)
plt.show()


