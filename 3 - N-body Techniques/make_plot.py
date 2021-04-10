from nbody import System
from sys import argv
import numpy as np
import matplotlib.pyplot as plt

N = len(argv) - 1
data = np.zeros((N, 8))
t = np.zeros(N)
for i in range(N):
	s = System.read(argv[i+1])
	
	t[i] = s.time
	sun = s.particles[0]
	for j in range(8):
		data[i, j] = sun.distance_from(s.particles[j+1])


for i in range(8):
	plt.plot(t, data[:,i])
	
plt.ylabel("radii (Au)")
plt.xlabel("time (years)")
plt.show()
	
	


