import numpy as np

#
# Uses Euler to solve the SHO ODE
#

def f(x, y, z):
	return z
	
def g(x, y, z):
	return -y

# step size
h = 0.01

# BCs
x0 = 0.0
y0 = 1.0
z0 = 0.0

# End point
xf = 50.0

# create arrays
x = np.arange(x0, xf, h)
N = len(x)

y = np.zeros(N)
y[0] = y0  # fixes BC
z = np.zeros(N)
z[0] = z0  # fixes BC

# now loop to advance the solution
for i in range(N-1):  # note the N-1 -- we're calculating the i+1th soltuion
	k = h * f(x[i], y[i], z[i])
	l = h * g(x[i], y[i], z[i])
	
	y[i+1] = y[i] + k
	z[i+1] = z[i] + l

# what's the actual solution?  Just cos(x):
y_actual = np.cos(x)

# print out the data so we can plot it later
np.savetxt("euler_data.txt", np.column_stack((x, y, y_actual)))

