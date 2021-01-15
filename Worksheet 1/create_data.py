import numpy as np

x = np.arange(0.0, 1.0, 0.01)
f = x*np.exp(-x*x)

np.savetxt("data.txt", np.column_stack((x, f)))
