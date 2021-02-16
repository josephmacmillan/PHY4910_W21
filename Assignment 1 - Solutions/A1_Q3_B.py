import numpy as np
import phy4910
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

# build a realistic model of a white dwarf
r, rho, R, M = phy4910.build_white_dwarf(3.789e6)

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

ax.plot(r, rho, label = f"{M:.3f} M$_\odot$")
ax.legend()
plt.savefig("A1_Q3_B.pdf")

