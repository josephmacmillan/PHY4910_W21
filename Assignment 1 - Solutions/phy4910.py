import numpy as np
    
"""
Solves a coupled pair of ODEs using euler
 
Takes as arguments:
  x_start - starting point for independent coordinate
  x_end - ending point
  h - difference between x_i and x_i+1 (i.e., delat x)
  y0 - initial value for first variable y(x_start)
  z_0 - initial value for second variable z(x_start)
  f - function for derivative of first variable (i.e., f = dy/dx)
  g - function for deriviative of second variable (i.e., g = dz/dx)
    
  returns three arrays, x[0,N-1], y[0,N-1], and z[0,N-1].
"""
def ode_euler(x_start, x_end, h, y0, z0, f, g):
    
    x = np.arange(x_start, x_end, h)
    N = len(x)
    y = np.zeros(N)
    y[0] = y0
    z = np.zeros(N)
    z[0] = z0
    
    for i in range(0, N-1):
        k1 = h * f(x[i], y[i], z[i])
        l1 = h * g(x[i], y[i], z[i])
            
        y[i+1] = y[i] + k1
        z[i+1] = z[i] + l1
                
    return x, y, z
    

""" 

Solves a coupled pair of ODEs using runge kutta.
 
Takes as arguments:
  x_start - starting point for independent coordinate
  x_end - ending point
  h - difference between x_i and x_i+1 (i.e., delat x)
  y0 - initial value for first variable y(x_start)
  z_0 - initial value for second variable z(x_start)
  f - function for derivative of first variable (i.e., f = dy/dx)
  g - function for deriviative of second variable (i.e., g = dz/dx)
  stop - a function for a stopping criteria.  Function must take three number (xi, yi, zi).
  data - any data that needs to be passed to the stop function in addition to xi, yi, zi.
  
  returns three arrays, x[0,N-1], y[0,N-1], and z[0,N-1].
"""

def ode_rk4(x_start, x_end, h, y0, z0, f, g, stop, data):
    
    x = np.arange(x_start, x_end, h)
    N = len(x)
    y= np.zeros(N)
    y[0] = y0
    z = np.zeros(N)
    z[0] = z0
    
    for i in range(0, N-1):
        k1 = h * f(x[i], y[i], z[i])
        l1 = h * g(x[i], y[i], z[i])
        k2 = h * f(x[i] + 0.5 * h, y[i] + 0.5 * k1, z[i] + 0.5 * l1)
        l2 = h * g(x[i] + 0.5 * h, y[i] + 0.5 * k1, z[i] + 0.5 * l1)
        k3 = h * f(x[i] + 0.5 * h, y[i] + 0.5 * k2, z[i] + 0.5 * l2)
        l3 = h * g(x[i] + 0.5 * h, y[i] + 0.5 * k2, z[i] + 0.5 * l2)
        k4 = h * f(x[i] + h, y[i] + k3, z[i] + l3)
        l4 = h * g(x[i] + h, y[i] + k3, z[i] + l3)
    
        y[i+1] = y[i] + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        z[i+1] = z[i] + (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0
        
        if stop(x[i+1], y[i+1], z[i+1], data):
            return x[0:i], y[0:i], z[0:i]
        
    return x, y, z
    
""" 
Builds a polytropic model white dwarf.
 
Takes as arguments:
  n - index of polytrope
  k - pressure-density constant in cgs units
  rho_c - central density in cgs units
  
  returns four arrays, eta (the dimensionless radii), varrho (the dimensionless density), 
          r (the physical radii in km), and rho (the physical density in g/cc)
"""
def build_white_dwarf_polytrope(n, k, rho_c):
	
	G = 6.6743e-8  # cgs units
	M_sun = 1.989e33 # grams
	R_sun = 6.957e10 # cm
	
	def f(eta, varrho, sigma):
		return sigma
		
	def g(eta, varrho, sigma):
		return -2.0/eta * sigma - np.power(varrho, n)
		
	def stop(eta, varrho, sigma, error):
		if abs(varrho) < error:
			print(f"### Stopping ODE solving at radius eta = {eta}")
			return True
		else:
			return False
		
	# solve the Lane-Emden equation using Runge-Kutta
	eta, varrho, sigma = ode_rk4(0.000000001, 15.0, 0.001, 1.0, 0.0, f, g, stop, 1e-7)

	# where is the surface?  We can find it by keeping only the positive values of varrho.
	# Warning: this assumes that the density function doesn't later on go positive again (i.e. oscillate)
	condition = varrho > 0.0
	eta = eta[condition]
	varrho = varrho[condition]

	surface = eta[-1]
	print(f"You built a polytropic white dwarf with index n = {n}")
	print(f"  * The dimensionless surface is at eta = {surface:.3f}")
	# what is the physical radius?
	lam = np.sqrt( (n+1) * k * np.power(rho_c, (1.0 - n) / n) / 4.0 / np.pi / G )  # cgs
	lam_km = lam/100000
	R = surface * lam_km
	print(f"  * The physical surface is at r = {R:.3f} km")

	# now do the integration
	m = np.trapz(np.power(varrho, n) * eta**2, eta)
	print(f"  * The dimensionless mass is {m:.3f}")

	# what is the physical mass of the white dwarf?
	M = 4.0 * np.pi * rho_c * lam**3 * m  #g
	M_solar = M / M_sun
	print(f"  * The physical mass is {M_solar:.3f} M_sun")


	# now print out the actual density data
	r = lam_km * eta
	rho = rho_c * np.power(varrho, n)
	
	return eta, varrho, r, rho

""" 
Builds a more realistic model of a white dwarf.
 
Takes as arguments:
  
  rho_c - central density in cgs units
  h - the step-size, in dimensionless radii (eta)
  eta_f - the max size of the eta array
  
  returns two arrays, r (the physical radii in km), and rho (the physical density in g/cc)
          and two numbers (radius of surface in km, and mass of white dwarf in solar masses)
"""
def build_white_dwarf(rho_c, h = 1e-3, eta_f = 20.0):

	# some useful constants
	G = 6.6743e-8  # cgs units
	M_sun = 1.989e33 # grams
	R_sun = 6.957e10 # cm
	
	k_nr = 3.166e12  # cgs; see Worksheet 3
	rho_0 = 3.789e6 # g/cc
	lam = np.sqrt( k_nr / (4.0 * np.pi * G * np.power(rho_0, 1.0/3.0)) ) # cm
	
	# initial condition is varrho_c; it's rho_c / rho_0.  See Assignment for details.
	varrho_c = rho_c / rho_0
	
	def A(x):
		return -(5.0/9.0) * np.power(x, -4.0/3.0) * np.power(1.0 + np.power(x, 2.0/3.0), -0.5) - (2.0/3.0) * np.power(x, -2.0/3.0) * np.power(1.0 + np.power(x, 2.0/3.0), -1.5) + (1.0/3.0) * np.power(1.0 + np.power(x, 2.0/3.0), -2.5)

	def B(x):
		return (5.0/3.0) * np.power(x, -1.0/3.0) * np.power(1.0 + np.power(x, 2.0/3.0), -0.5) - (1.0/3.0) * np.power(x, 1.0/3.0) * np.power(1.0 + np.power(x, 2.0/3.0), -1.5)
        
	def f(eta, varrho, sigma):
		return sigma

	def g(eta, varrho, sigma):
		return -2.0/eta * sigma - A(varrho) / B(varrho) * sigma*sigma - 1.0/B(varrho) * varrho

	def stop(eta, varrho, sigma, data):
		arrho_c = data[0]
		rel_error = data[1]
		if varrho/varrho_c < rel_error:
			print(f"### Stopping ODE solving at radius eta = {eta}")
			return True
		else:
			return False

	# solve the ODE using Runge-Kutta
	print(f"Building a white dwarf; stand by ...")
	eta, varrho, sigma = ode_rk4(1e-5, eta_f, h, varrho_c, 0.0, f, g, stop, [varrho_c, 1e-6])
	print(f"  Surface at eta = {eta[-1]}")
	# calculate the dimensionless mass
	m = np.trapz(varrho * eta**2, eta)
	print(f"  Dimensionless mass = {m}")
    
	# this is all dimensionless stuff, so let's add the units in.
	# I prefer km to cm, though, for the radius, and 10^6 g/cm^3:
	
	r = eta * lam / 100000.0
	R = r[-1]
	rho = varrho * rho_0 / 1e6
	M = (4.0 * np.pi * rho_0 * lam**3 * m) / M_sun

	print(f"You built a white dwarf")
	print(f"  * Central density rho_c = {rho_c / 1e6}*10^6 g/cm^3")
	print(f"  * Radius R = {R:.0f} km")
	print(f"  * Mass M = {M:.3f} M_sun")

	return r, rho, R, M



