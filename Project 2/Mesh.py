# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 14:32:20 2021

@author: Darryen
"""
import numpy as np
import matplotlib.pyplot as plt
import nbody

class Mesh():
    '''
    A class for the mesh of our galaxy.
    '''
    
    def __init__(self, L, Ng, omegaM0 = 0.27, omegaL0 = 0.72, H0 = 71):
        self.L = L
        self.Ng = Ng
        self.OmegaM0 = omegaM0
        self.OmegaL0 = omegaL0
        self.H0 = H0
        
    def f(self, a):
        return (self.OmegaM0/a + (1 - self.OmegaM0 - self.OmegaL0) - self.OmegaL0*a**2)**(-1/2)
    
    def G(self, a, l, m, n):
        kx = (2*np.pi*l)/(self.Ng)
        ky = (2*np.pi*m)/(self.Ng)
        kz = (2*np.pi*n)/(self.Ng)
        Omega0 = self.OmegaL0 + self.OmegaM0
        
        if l == 0 and m == 0 and n == 0:
            return 0
        else:
            return (3*Omega0 /(8*a)) * (np.sin(kx/2)**2 + np.sin(ky/2)**2 + np.sin(kz/2)**2)**(-1)
    
    def calc_density(self, system):
        '''
        Parameters
        ----------
        system : System Object
    
        Returns
        -------
        Relativistic Density: Array of 3 
        '''
        print('Calculating the density of the system ... \n')
        
        # Start with zeroes
        rho = np.zeros((self.Ng, self.Ng, self.Ng))
        for particle in system.particles:
            # Takes the density in the xth, yth, zth direction
            i = int(particle.position[0] * self.Ng / self.L)
            j = int(particle.position[1] * self.Ng / self.L)
            k = int(particle.position[2] * self.Ng / self.L)
        
            # We need to consider boundary conditions with the following code
            if i == self.Ng:
                i = 0
            if j == self.Ng:
                j = 0
            if k == self.Ng:
                k = 0
        
            # Now we start to count the number of particles in the grid cell
            rho[i,j,k] += 1
            
            # What we really want is delta which is rho/rho_B0 - 1
            
            # First we need the number of particles along a direction
            N = len(system.particles)**(1/3)
            
            # Now if we assume the density is uniformly distributed 
            # (I think this is a fair assumption given the size of our simulation)
            # we can find the average number of particles in the cell
            N_avg = (N/self.Ng) ** 3
            
            # Now the relativistic density is
            delta = rho/N_avg - 1
            
            # and the average of that is 
            delta_avg = np.average(delta)
            
            # Let's print a few values for reference
            print(f"Average number of particles: {N_avg:.3f} \n  Average relativistic density: {delta_avg:.3f}")
            return delta
        
    
    def calc_potential(self, delta, a):
        '''
        Parameters
        ----------
        delta : Array of 3
                The density of a cell
        a : TYPE
            Scale factor

        Returns
        -------
        Potential: Array of 3
                   The potential of a cell 
        '''
        print("Calculating the potential from the density... \n")
        
        # Let's get the FFT of the density
        delta_t = np.fft.fftn(delta)
        
        # Now we need the fourier potential (this is going to be a complex array)
        phi_t = np.zeros((self.Ng, self.Ng, self.Ng), dtype = complex)
        
        # We are going to use a nested loop for each fourier direction (l,m,n)
        for l in range(self.Ng):
            for m in range(self.Ng):
                for n in range(self.Ng):
                    phi_t[l,m,n] = self.G(a,l,m,n)*delta_t[l,m,n]
                    
        # An inverse fourier transform to get the potential we need
        phi = np.fft.ifftn(phi_t)
        
        # Take only the real component of the potential 
        return phi.real
    
    def calc_accels(self, system, a):
        '''
        Parameters
        ----------
        system : System object
            A system of particles
        a : Float
            The scale factor
        Returns
        -------
        None.
        
        Accelerates and updates the particles
        '''
        # We need delta and phi
        delta = self.calc_density(system)
        phi = self.calc_potential(delta, a)
        
        # Initialize each of the accelerations to zero
        gx = np.zeros((self.Ng, self.Ng, self.Ng))
        gy = np.zeros((self.Ng, self.Ng, self.Ng))
        gz = np.zeros((self.Ng, self.Ng, self.Ng))
        
        # Now we update the accelerations
        for i in range(self.Ng - 1):
            for j in range(self.Ng - 1):
                for k in range(self.Ng - 1):
                    gx[i,j,k] = (-1/2)*(phi[i+1,j,k] - phi[i-1,j,k])
                    gy[i,j,k] = (-1/2)*(phi[i,j+1,k] - phi[i,j-1,k])
                    gz[i,j,k] = (-1/2)*(phi[i,j,k+1] - phi[i,j,k-1])
        
        
        
        # Extra code that puts each direction of acceleration into 1 array
        accel = np.zeros((len(system.particles),3))
        
        for p in range(len(system.particles)):
            # Takes each direction 
            i = int((system.particles[p].position[0] * self.Ng / self.L))
            j = int((system.particles[p].position[1] * self.Ng / self.L))
            k = int((system.particles[p].position[2] * self.Ng / self.L))
            
            if i == self.Ng:
                i = 0
            if j == self.Ng:
                j = 0
            if k == self.Ng:
                k = 0
            
            accel[p] = [gx[i,j,k], gy[i,j,k], gz[i,j,k]]
            
            return accel
    
    def move_particles(self, system, accels, a, da):
        '''
        Parameters
        ----------
        system : System object
                A system input
        a : Float
            Scale Factor
        da : Float
            Change in scale factor
        Returns
        -------
        Updates the particles positions and velocities
        '''
        # Find all the particles in a system
        for i in range(len(system.particles)):
            # Finds the particles across 3 dimensions
            for j in range(3):
                system.particles[i].velocity[j] += self.f(a) * accels[i][j] * da
                system.particles[i].position[j] += self.f(a)/a**2 * system.particles[i].velocity[j]*da
                
                # Checks to make sure the periodic boundary conditions are satisfied
                while system.particles[i].position[j].any() > self.L:
                    system.particles[i].position[j] -= self.L
                    
                while system.particles[i].position[j].any() < 0.0:
                    system.particles[i].position[j] += self.L
                    
                    
def main():
    return 0
    
        
if __name__ == '__main__':
    main()
    
        