import numpy as np
import matplotlib.pyplot as plt

z = np.loadtxt('data.txt', unpack=True)

plt.plot(z[0],z[1])
plt.show()