import numpy as np
from PartA import *
import scipy.integrate as integ
from sys import argv
from numpy.fft import fftn, ifftn
from nbody import *



class Mesh:
    def __init__(self, L, Ng, OmegaM0 = OmeM0, OmegaL0 = OmeL0, H0 = H_0, rho_c = rho_c, rho_m = rho_m0, rho_L = rho_L):
        
        self.L = L
        
        self.Ng = Ng
        self.H0 = H0
        self.OmegaL0 = OmegaL0
        self.OmegaM0 = OmegaM0
        self.OmegaCu = 1 - OmegaM0 - OmegaL0
        
        self.rho_c = rho_c
        self.rho_m = rho_m
        self.rho_L = rho_L
        
        print(f"Length = {self.L}")
        print(f"Number of cells = {self.Ng}")
        print(f"H0 = {self.H0}")
        print(f"OmegaM0 = {self.OmegaM0}")
        print(f"OmegaL0 = {self.OmegaL0}")
        print(f"rho_c = {self.rho_c}")
        print(f"rho_m = {self.rho_m}")
        print(f"rho_L = {self.rho_L}")
        
    def f(self,a):
        f = 1/np.sqrt((self.OmegaM0/a) + (1 - self.OmegaM0 - self.OmegaL0) - (self.OmegaM0**2))
        #print(f"f = {f}")
        return f
        
    def g(self,a):
        g = self.f(a)/(a**2)
        #print(f"g = {g}")
        return g
    
    def G(self, l, m, n, a):
        k_0 = 2*np.pi/self.Ng
        if l==0 and m==0 and n==0:
            G = 0
        else:
            G = -3*self.OmegaM0 / (8*a* ((np.sin(k_0*l/2)**2) + (np.sin(k_0*m/2)**2) + (np.sin(k_0*n/2)**2)))
        #print(f"G(k) = {G}")
        return G
        
    def calc_density(self, parts,a):
        rho = np.zeros((self.Ng, self.Ng, self.Ng))
        
        for part in parts.particles:
            i = int(part.position[0] * self.Ng / self.L)
            j = int(part.position[1] * self.Ng / self.L)
            k = int(part.position[2] * self.Ng / self.L)
            
            rho[i,j,k] += 1
        
        Nt = np.power(len(parts.particles),1/3)
        #print(f"Total Number of Particles = {Nt}")
        #print(f"No. of cells = {self.Ng}")
        Navg = np.power(Nt/(self.Ng),3)
        #print(f"Average number of Particles = {Navg}")
        delta = rho/Navg - 1
        
        return delta
    
    def calc_potential(self, parts, a):
        
        delta = self.calc_density(parts, a)
        delta_t = fftn(delta)
        #print(f"Shape of Delta_t = {np.shape(delta_t)}")
        phi_t = np.zeros((self.Ng,self.Ng,self.Ng), dtype=complex)
        #print(f"Shape of phi_t = {np.shape(phi_t)}")
        for l in range(self.Ng):
            for m in range(self.Ng):
                for n in range(self.Ng):
                    g = self.G(l, m, n, a)
                    phi_t[l,m,n] = g * delta_t[l,m,n]
                    
        phi_unreal = ifftn(phi_t)
        phi = phi_unreal.real
        return phi
        
    
    def calc_accels(self,parts,a):
        delta = self.calc_density(parts, a)
        phi = self.calc_potential(parts, a)
        #print(f"Shape of phi = {np.shape(phi)}")
        
        g_x = np.zeros((self.Ng,self.Ng,self.Ng))
        g_y = np.zeros((self.Ng,self.Ng,self.Ng))
        g_z = np.zeros((self.Ng,self.Ng,self.Ng))
        #print(f"THIS g_accel SIZE = {np.shape(g_x)}")
        
        for i in range(self.Ng):
            for j in range(self.Ng):
                for k in range(self.Ng):
                    
                    g_x[i,j,k] = -0.5*(phi[(i+1)%self.Ng,j,k] - phi[(i-1)%self.Ng,j,k])
                    g_y[i,j,k] = -0.5*(phi[i,(j+1)%self.Ng,k] - phi[i,(j-1)%self.Ng,k])
                    g_z[i,j,k] = -0.5*(phi[i,j,(k+1)%self.Ng] - phi[i,j,(k-1)%self.Ng])
                    
                    
        accel = np.zeros((len(parts.particles),3))
        #print(f"Shape of accelarr = {np.shape(accel)}")
        #print(f"Shape of g = {np.shape(g_x)}")
        for n in range(len(parts.particles)):
            i = int((parts.particles[n].position[0] * self.Ng / self.L))
            j = int((parts.particles[n].position[1] * self.Ng / self.L))
            k = int((parts.particles[n].position[2] * self.Ng / self.L))
            
            if i == self.Ng:
                i = 0
            if j == self.Ng:
                j = 0
            if k == self.Ng:
                k = 0
            
            accel[n] = [g_x[i,j,k], g_y[i,j,k], g_z[i,j,k]]
        print(f"Acel TEST = {accel[0]}")
        #print(f"THIS ACCEL SIZE = {np.shape(accel)}")
        #print(f"Accleration Array = {accel}")
        
        return accel
    
    def move_particles(self,parts, a, da, accel):
        
        for i in range(len(parts.particles)):            
            for j in range(3):
                #print(f"Iteration = {i}")
                #print(f"Size of accels = {np.shape(accel)}")
                #print(f"Particle Velocity = {parts.particles[i].velocity[j]}")
                parts.particles[i].velocity[j] += self.f(a) * accel[i][j] * da
                parts.particles[i].position[j] += self.g(a) * parts.particles[i].velocity[j] * da 
                while parts.particles[i].position[j] > self.L:
                    parts.particles[i].position[j] -= self.L
                while parts.particles[i].position[j] < 0:
                    parts.particles[i].position[j] += self.L
                
            
            
def update(a,da):
    """
    

    Parameters
    ----------
    a : Scale factor - Start 'time'
    da : 'Timestep' - Scale factor interval size

    Returns
    -------
    Writes updated Particle attributes to data files over scale factor ('time')

    """
    s = System.read_binary("init_z30_L50_N64.dat")
    
    m = Mesh(50,128)
    count = 1
    while a <= 1:
        accels = m.calc_accels(s,a)
        
        m.move_particles(s,a,da,accels)
        
        s.write(f"./ori/Naru/Naru_{a:5.3f}.dat")
        print(f"File {count} created")
        count +=1
        a += da

a = 1/31
n = 100
da = (1 - a)/n        
update(a,da)