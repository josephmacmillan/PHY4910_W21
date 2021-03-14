import numpy as np
import sys

#creates a set of 20 numbers from -1 to 1500 of log base 10

m = np.logspace(-1, np.log10(1500), 20)
if len(sys.argv) > 1:
    i = int(sys.argv[1])
    if i < 20 and i >= 0:
        print(m[i])
        sys.exit()

print(m)

