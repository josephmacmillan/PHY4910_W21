from math import sqrt, pi
from nbody import Particle

def calc_energy(p1, p2):
	
	T = 0.5 * p1.mass * p1.velocity2() + 0.5 * p2.mass * p2.velocity2()

	r = p1.distance_from(p2)
	U = -p1.mass * p2.mass / r
	
	E = T + U 

	return E

def calc_momentum(p1, p2):
	
	px = p1.mass * p1.vx() + p2.mass * p2.vx()
	py = p1.mass * p1.vy() + p2.mass * p2.vy()
	pz = p1.mass * p1.vz() + p2.mass * p2.vz()
	
	return sqrt(px**2 + py**2 + pz**2)
	

def evolve(p1, p2, dt, tf, filename):

	f = open(filename, "w")
	
	E0 = calc_energy(p1, p2)
	P0 = calc_momentum(p1, p2)
	
	print(f"# Evolving two-body system; Ei = {E0}, Pi = {P0}")
	
	t = 0
	while t < tf:
	
		r = p1.distance_from(p2)
		r3 = r*r*r
	
		for i in range(3):
			
			p1.accel[i] = -p2.mass * (p1.position[i] - p2.position[i]) / r3
			p2.accel[i] = -p1.mass * (p2.position[i] - p1.position[i]) / r3
			
			p1.velocity[i] += p1.accel[i] * dt
			p2.velocity[i] += p2.accel[i] * dt
			
			p1.position[i] += p1.velocity[i] * dt
			p2.position[i] += p2.velocity[i] * dt
			
		t += dt
		
		E = calc_energy(p1, p2)
		P = calc_momentum(p1, p2)
		
		f.write(f"{t} {p1.x()} {p1.y()} {p1.z()} {p2.x()} {p2.y()} {p2.z()}\n")
	
	f.close()
	
	print(f"Done! Ef = {E}, Pf = {P}")


earth = Particle(0.1, [0.9, 0, 0], [0, -1.2, 0])
sun = Particle(0.9, [-0.1, 0, 0], [0, 0.1333333333333333, 0])
print(calc_energy(earth, sun))
print(calc_momentum(earth, sun))

evolve(earth, sun, 0.01, 500, "earthsun.txt")

