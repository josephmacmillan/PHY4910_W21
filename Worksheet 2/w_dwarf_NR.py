import numpy as np
import matplotlib.pyplot as plt
import phy4910_rk4 as rk4


#ODE Functions
"""
x = eta, radius (dimensionless?)
y = rho, density (definitely dimensionless)
z = sigma, uhh ??
"""
n = 1.5
def f(x,y,z):
    return z

def g(x,y,z):
    func = (-(2*z)/x)-(y**n)
    return func


#Initial values
x_start = 0.0000001  #zero returns error, so this is close to zero
x_stop = 4

h = 0.0001  #steps for rk4 ODEsolver
y0 = 1
z0 = 0


#constants and defined values
k = 3.166 * (10**12)    #polytrope constant
pc = 4.045 * (10**6)    #rho_c, denstiy at center... centre?
G = 6.6743 * (10**(-8)) #cgs


#ODE Solver
x,y,z = rk4.ode_solver(x_start, x_stop, h, y0, z0, f, g)
condition = y > 0.0     #Ignores negative y values
x = x[condition]
y = y[condition]


#Mass integral
m = (y**n)*(x**2)
mass = np.trapz(m,x)

#Radial Scale Factor
lmda = np.sqrt(((n+1)*k*(pc**( (1-n)/n)))/((4*np.pi*G))) 
x_max = lmda*x[-1]*(10**-5)  #last value of x converted cm to km [and not m]


#Density and Radius with dimensions 
y_cgs = y*pc        #g/cm^3
x_km = x*lmda*(10**-5)       #km

#Mass with dimensions
M = 4*np.pi*pc*(lmda**3)*mass*(10**-3)   #AAHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
SM = M/(2*(10**30))          #converted to Solar Mass


#Plots 
plt.plot(x,y)
plt.xlabel('Density')
plt.ylabel('Radius')
plt.grid(True, axis="both")
plt.show()

plt.plot(x_km, y_cgs)
plt.xlabel('Radius (km)')
plt.ylabel('Density (g/cm$^3$)')
plt.grid(True, axis="both")
plt.show()

print('surface radius=', x[-1])
print('mass=', mass)
print('lambda=',lmda)
print('radius=',x_max,'km')
print('Mass =',M,'kg, or', SM, 'M\u2609')