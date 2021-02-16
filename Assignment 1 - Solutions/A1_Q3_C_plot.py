import numpy as np
import phy4910
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':14})
plt.rc('axes', labelsize=16)
plt.rcParams.update({'figure.autolayout': True})

(rho_c, R, M) = np.loadtxt("A1_Q3_C_data.txt", unpack=True)

# first plot - i. - M versus log(rho_c)
fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
ax.grid(True)
ax.set_xlabel(r"$\rho_c$ ($10^6$ g/cm$^3$)")
ax.set_ylabel(r"$M$ (M$_\odot$)")
ax.semilogx(rho_c/1e6, M)
plt.savefig("A1_Q3_C_i.pdf")

# second plot - i. - R versus log(rho_c)
fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
ax.grid(True)
ax.set_xlabel(r"$\rho_c$ ($10^6$ g/cm$^3$)")
ax.set_ylabel(r"$R$ (km)")
ax.semilogx(rho_c/1e6, R)
plt.savefig("A1_Q3_C_ii.pdf")

# third plot - i. - M versus R
fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(1,1,1)
ax.grid(True)
ax.set_ylabel(r"$M$ (M$_\odot$)")
ax.set_xlabel(r"$R$ (km)")
ax.set_xlim(0, 25000)
ax.set_ylim(0, 1.5)
ax.plot(R, M, color="black")
ax.plot(R, M[0] * R[0]**3 / (R**3), "--", color="black")
ax.plot(R, np.full(len(rho_c), M[-1]),  "--", color="black")
plt.savefig("A1_Q3_C_iii.pdf")
