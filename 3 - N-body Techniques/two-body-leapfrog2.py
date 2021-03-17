from math import sqrt, pi
from nbody import Particle, System

earth = Particle(0.1, [0.9, 0.0, 0.0], [0.0, -1.2, 0.0])
sun = Particle(0.9, [-earth.mass / 0.9 * earth.x(), 0.0, 0.0], [0.0, -earth.mass / 0.9 * earth.vy(), 0.0])

s = System([sun, earth])

dt = 0.01
tf = 500
ti = 0

f = open("two-body-leap2.dat", "w")

t = ti
while t < tf:

	s.update_positions(0.5 * dt)
	
	s.calc_accels()
	s.update_velocities(dt)
	s.update_positions(0.5 * dt)
	
	t += dt
	s.time = t

	#print(f"t = {t} | dE = {(s.T + s.U - E0)/E0}")
	f.write(f"{t} {s.particles[0].x()} {s.particles[0].y()} {s.particles[0].z()} {s.particles[1].x()} {s.particles[1].y()} {s.particles[1].z()}\n")

f.close()

print(f"# Done! Ef = {s.T + s.U}")

