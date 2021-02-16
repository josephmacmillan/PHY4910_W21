import numpy as np
import phy4910

N = 25
rho_c = np.logspace(4, 12, N)
M = np.zeros(N)
R = np.zeros(N)

for i in range(N):
	r, rho, R[i], M[i] = phy4910.build_white_dwarf(rho_c[i])
	
np.savetxt("A1_Q3_C_data.txt", np.column_stack((rho_c, R, M)))

