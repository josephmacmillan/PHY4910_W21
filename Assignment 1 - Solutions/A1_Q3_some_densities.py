import numpy as np
import phy4910
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

# build some realistic models of white dwarfs and plot the densities

rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':14})
plt.rc('axes', labelsize=16)
plt.rcParams.update({'figure.autolayout': True})

fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
ax.grid(True)
#ax.set_xlim(8e-5, 4e-2)
#ax.set_ylim(8e-4, 2e4)
ax.set_xlabel(r"$r$ (km)")
ax.set_ylabel(r"$\rho$ ($10^6$ g/cm$^3$)")

rho_0s = [3.789e4, 3.789e6, 3.789e8, 3.789e10, 3.789e12]
colors = ['black', 'blue', 'red', 'green', 'orange']

for i in range(len(rho_0s)):
	r, rho, R, M = phy4910.build_white_dwarf(rho_0s[i], h=1e-4)
	ax.loglog(r, rho, color=colors[i], label=f"{M:.3f} M$_\odot$")

ax.legend()
#plt.show()
plt.savefig("A1_Q3_some_densities.pdf")
	
	
