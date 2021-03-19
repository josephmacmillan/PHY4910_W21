from math import sqrt, pi
from nbody import Particle, System


s = System.read("solar_system.dat")

# calculate energies for later
s.calc_accels()
E0 = s.T + s.U

# units are now in AUs, years, and earth masses
dt = 0.001
tf = 200.0
ti = 0

dt_out = 1
dt_log = 0.1

t = ti
t_out = ti
t_log = ti
while t < tf:

	s.update_positions(0.5 * dt)
	
	s.calc_accels()
	s.update_velocities(dt)
	s.update_positions(0.5 * dt)
	
	t += dt
		
	s.time = t
	
	if t >= t_out:
		fname = f"run_{t:07.2f}.dat"
		s.write(fname)
		t_out += dt_out
        
	if t >= t_log:
		print(f"t = {t} | dE = { (s.T + s.U - E0)/E0}")
		t_log += dt_log

print(f"Done!")

