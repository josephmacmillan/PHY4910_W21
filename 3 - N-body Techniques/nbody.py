from math import sqrt

class Particle:

	def __init__(self, mass, pos, vel):
	
		self.mass = mass
		self.position = pos
		self.velocity = vel
		
		self.accel = [0.0, 0.0, 0.0]
		
	def __str__(self):
		return f"{self.mass} {self.position[0]} {self.position[1]} {self.position[2]} {self.velocity[0]} {self.velocity[1]} {self.velocity[2]}"
		
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
		dx = p.x() - self.x()
		dy = p.y() - self.y()
		dz = p.z() - self.z()
		return sqrt(dx**2 + dy**2 + dz**2)

	def abs_velocity(self):
		return sqrt(self.velocity[0]**2 + self.velocity[1]**2 + self.velocity[2]**2)

	def velocity2(self):
		return (self.velocity[0]**2 + self.velocity[1]**2 + self.velocity[2]**2)
		
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






