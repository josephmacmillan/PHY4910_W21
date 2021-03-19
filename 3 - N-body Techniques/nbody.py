from math import sqrt

class Particle:

	def __init__(self, mass, pos, vel):
		self.mass = mass
		self.position = pos
		self.velocity = vel

		self.accel = [0.0, 0.0, 0.0]

	def x(self):
		return self.position[0]

	def y(self):
		return self.position[1]

	def z(self):
		return self.position[2]

	def vx(self):
		return self.velocity[0]

	def vy(self):   
		return self.velocity[1]

	def vz(self):
		return self.velocity[2]

	def radius(self):
		return sqrt(self.position[0]**2 + self.position[1]**2 + self.position[2]**2)

	def distance_from(self, p):
		return sqrt( (self.position[0] - p.position[0])**2 + (self.position[1] - p.position[1])**2 + (self.position[2] - p.position[2])**2 )

	def abs_velocity(self):
		return sqrt(self.velocity[0]**2 + self.velocity[1]**2 + self.velocity[2]**2)

	def velocity2(self):
		return (self.velocity[0]**2 + self.velocity[1]**2 + self.velocity[2]**2)

	def __str__(self):
		return f"{self.mass} {self.position[0]} {self.position[1]} {self.position[2]} {self.velocity[0]} {self.velocity[1]} {self.velocity[2]}"

	@classmethod
	def from_string(cls, s):
		parts = s.split()
		m = float(parts[0])
		x = float(parts[1])
		y = float(parts[2])
		z = float(parts[3])
		vx = float(parts[4])
		vy = float(parts[5])
		vz = float(parts[6])

		return cls(m, [x,y,z], [vx,vy,vz])

class System:

	def __init__(self, parts, t = 0.0, G = 1.0):
		self.N = len(parts)
		self.time = t
		self.particles = parts

		self.T = 0.0
		self.U = 0.0

		self.G = G

	def __str__(self):
		ret = f"{self.N}\n{self.time}\n{self.G}\n"
		for i in range(self.N):
			ret += str(self.particles[i]) + '\n'
		return ret

	def write(self, filename):
		f = open(filename, 'w')
		print(f"Writing to file {filename} ({self.N} particles at time {self.time})")
		f.write(str(self))
		f.close()

	@classmethod
	def read(cls, filename):
		f = open(filename, 'r')
		N = int(f.readline())
		t = float(f.readline())
		G = float(f.readline())
		print(f"Reading from file {filename} ({N} particles at time {t}")
		parts = []
		for i in range(N):
			parts.append(Particle.from_string(f.readline()))
		f.close()
		return cls(parts, t, G)    

	def all_x(self):
		x = []
		for i in range(self.N):
			x.append(self.particles[i].x())

		return x

	def all_y(self):
		y = []
		for i in range(self.N):
			y.append(self.particles[i].y())

		return y    

	def all_z(self):
		z = []
		for i in range(self.N):
			z.append(self.particles[i].z())

		return z

	def calc_accels(self):
	
		for i in range(self.N):
			self.particles[i].accel = [0,0,0]
			self.U = 0.0
			self.T = 0.0

		for i in range(self.N):
			p_i = self.particles[i]
			for j in range(i+1, self.N):
				p_j = self.particles[j]
				r = p_i.distance_from(p_j)
				r3 = r**3
				for k in range(3):
					p_i.accel[k] += -self.G * p_j.mass * (p_i.position[k] - p_j.position[k]) / r3
					p_j.accel[k] += -self.G * p_i.mass * (p_j.position[k] - p_i.position[k]) / r3
				self.U -= self.G * p_j.mass * p_i.mass / r
			self.T += 0.5 * p_i.mass * p_i.velocity2()

	def update_positions(self, dt):
		for i in range(self.N):
			for k in range(3):
				self.particles[i].position[k] += self.particles[i].velocity[k] * dt

	def update_velocities(self, dt):
		for i in range(self.N):
			for k in range(3):
				self.particles[i].velocity[k] +=  self.particles[i].accel[k] * dt



if __name__ == "__main__": 

	#p1 = Particle(1.0, [2, 3, 4], [5, 6, 7])
	#p2 = Particle.from_string("8 9 10 11 12 13 14")

	#s = System([p1, p2], t = 1231245)
	#print(s)
	#s.write("test.dat")
	
	s = System.read("test.dat")
	print(s)
