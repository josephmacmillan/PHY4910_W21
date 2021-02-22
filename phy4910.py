"""
Tech of modern astro
team Planck
Worksheet 2

Author: Georges Karagozian
"""

import numpy as np

#Solves ODE in terms of x (eta), y (rho), and z (sigma) values
def solveODE(x_start, x_stop, h, y0, z0, f, g, Type = "rk"):
   
    print(f"ODE Solver is Running at {h} increments.")
    #initialize the arrays for each variable
    x = np.arange(x_start,x_stop,h)
    N = len(x)
    y = np.zeros(N)
    z = np.zeros(N)
    
    #set initial conditions
    y[0] = y0
    z[0] = z0
    
    # if conditional statement to decide which solution method to use
    
    # default to runge-kutta
    if Type == "rk" :
        for i in range(0,N-1):
            # Loops to find new y and z based on functions and ki/li values
            k1 = h * f(x[i],y[i],z[i])
            l1 = h * g(x[i],y[i],z[i])
            
            k2 = h * f(x[i] + h/2, y[i] + k1/2, z[i] + l1/2)
            l2 = h * g(x[i] + h/2, y[i] + k1/2, z[i] + l1/2)
            
            k3 = h * f(x[i] + h/2, y[i] + k2/2, z[i] + l2/2)
            l3 = h * g(x[i] + h/2, y[i] + k2/2, z[i] + l2/2)
            
            k4 = h * f(x[i] +h, y[i] +k3, z[i] +l3)
            l4 = h * g(x[i] +h, y[i] +k3, z[i] +l3)
            
            
            #use kn's and ln's (derivatives) to compute yn+1 and zn+1 (following points)
            y[i+1] = y[i] + k1/6 + k2/3 + k3/3 + k4/6
            z[i+1] = z[i] + l1/6 + l2/3 + l3/3 + l4/6
    
    # euler method
    elif Type == "eu" :
        
        for i in range(0, N-1):
            
            y[i+1] = y[i] + z[i] * h
            z[i+1] = z[i] - y[i] * h
    
    # in case a random string is passed display error    
    elif Type != "eu" and Type != "rk":
        print("Error! Please enter a valid type of solution method.")
    
    # return arrays of x,y,z variables
    return x,y,z    
        
        
#plots        
def plot(x,y, xlabel, ylabel, name = " "):
    import matplotlib.pyplot as plt
    
    #plot the arrays
    plt.plot(x,y)
    

    #label axes and make a title for the plot
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    #save plot if indicated in arguments
    if name != " ":
        plt.savefig(name)
    #display the plot
    plt.show()
        
    return None    
        
 

#Behavior or non-relativistic white dwarf star
def Nonrel_WhiteDwarf(x_stop, rho_c):
    #import libraries
    import numpy as np
    
    # define global variables for our constants
    n = 1.5
    k = 3.166 * 10**(12)
    G = 6.6743 * 10**(-8)
    #define functions in our problem in terms of eta,rho, and sigma
    #eta= the radial distance, rho= the density, sigma= derivative of rho
    def f(eta,rho,sigma):
        
        return sigma
    
    def g(eta,rho,sigma):
        
        return ((-2*sigma)/eta) - rho**n
    #call solveODE to solve for values of eta, rho, sigma
    eta,rho,sigma = solveODE(0.000001,x_stop,0.001,1,0,f,g)

    #trim the part of the arrays that's not a number
    condition = rho > 0.0
    eta = eta[condition]
    rho = rho[condition]
    sigma = sigma[condition]
    
    
   
    
    #Solving lambda, convert SI units
    lamb = (((n + 1) * k * rho_c**((1-n)/n))/ (4 * np.pi * G))**(0.5) #in cm
    lamb = lamb * 10**(-5)  #convert to km
    
    #calculate dimensionless mass
    y = rho**n * eta**2
    m = np.trapz(y,eta)
    
    #Calculate mass M, and convert to SI units (kg)
    M = 4 * np.pi * rho_c * (lamb * 10**(5))**3 * m * 10**(-3)
      
    
    #find real density and radius of the star
    eta = eta * lamb  # convert to km using the constant lambda
    rho = rho * rho_c 
    radius = eta[-1]
     
    print(f"White Dwarf Mass: {M}, Radius: {radius}, Density: {rho[0]}")
        
    return eta, rho, M, radius


def rel_WhiteDwarf(x_stop, rho_c):
    #import libraries
    import numpy as np
    
    # define variables for our constants
    n = 3
    k_r = 2.936 * 10**(14)
    G = 6.6743 * 10**(-8)
    
    def f(eta,rho,sigma):
    
        return sigma

    def g(eta,rho,sigma):
        
        return ((-2*sigma)/eta) - rho**n
    
    #call solveODE to solve for values of eta, rho, sigma
    eta,rho,sigma = solveODE(0.00001,x_stop,0.01,1,0,f,g)
    
    #trim the part of the arrays that's not a number
    condition = rho > 0.0
    eta = eta[condition]
    rho = rho[condition]
    sigma = sigma[condition]
    
    #calculate dimensionless mass
    y = rho**n * eta**2
    m = np.trapz(y,eta)

    #write formula for lambda then convert to km
    lamb = (((n + 1) * k_r * rho_c**((1-n)/n))/ (4 * np.pi * G))**(0.5) #in cm
    lamb = lamb * 10**(-5)  #convert to km
        
    #find real density and radius of the star

    eta = eta * lamb  # km
    rho = rho * rho_c
    
    #Calculate mass M
    #within the formula multiply lambda by 10^5 to convert to cm
    # our final asnwer for M will be in grams, so we will multiply the entire 
    #formula by 10^-3 to convert to kg
    
    M = 4 * np.pi * rho_c * (lamb * 10**(5))**3 * m * 10**(-3)
   
    
    radius = eta[-1]
    print(f"White Dwarf Mass: {M}, Radius: {radius}, Density: {rho[0]}")
    return eta, rho, M, radius




#bin a set of data Return horizontal axis (min to max of data) and binned data
def bin_data(data, nbins):
    
    #create an array filled with zeros with length nbins
    bins = np.zeros(nbins)
    
    dmax = np.amax(data)
    dmin = np.amin(data)
    
    #Sorts values from data into bins
    #b determines which bin the current number goes in and adds 1 to the bin it falls into (see lecture 3 for b equation)
    for i in range(len(data)):
    
        b = int((data[i]-dmin)/(dmax + 1e-10 - dmin)*nbins)
        bins[b] += 1
        
    #an array that sorts data from dmin to dmax with length nbins    
    horizontal_axis = np.linspace(dmin,dmax,nbins)
    print(f"The min/max data is {horizontal_axis}")
        
    return horizontal_axis, bins




#return random theta and phi values 
rng = np.random.default_rng()

def pick_direction():
    
    #correction for the distribution to prevent clustering at the peaks 
    #see link: http://corysimon.github.io/articles/uniformdistn-on-sphere/
    theta = np.arccos(1 - 2 * rng.random())

    theta = np.pi*rng.random() 
    phi = 2*np.pi*rng.random() 
    
    print(f"your random numbers are; {theta}, and {phi}")
    return theta, phi




#return random optical depth distance (tau)
def pick_optical_depth():
    x = rng.random()
    tau = -np.log(1-x)
    
    return tau
    print(f"Random optical depth distance is {tau}")


#models the path of a photon moving though an atmosphere. Returns two angles and the last scattering direction
def move_photon(tmax, zmax):
    x = []
    y = []
    z = []
    
    #starting at the origin
    x.append(0)
    y.append(0)
    z.append(0)
    
    #loop until photon moves from origin to the surface
    while True: 
        theta, phi = pick_direction()
        opdepth = pick_optical_depth()
        
        
        #how far the photon moves in the atmophere
        s = opdepth/tmax
        
        #how far the photon has moved since last time (not current position in atmosphere)
        dx = s * np.sin(theta)*np.cos(phi)
        dy = s * np.sin(theta)*np.sin(phi)
        dz = s * np.cos(theta)
    
        #last scattering direction of the photon at the surface
        x.append(x[-1] + dx)
        y.append(y[-1] + dy)
        z.append(z[-1] + dz)
        
        #Tell user the information
        print(f"Photon Position{x[-1]},{y[-1]},{z[-1]}")
        
        #if photon happens to go in the wrong direction, reset position of the photon
        if z[-1] < 0:
            x.clear()
            y.clear()
            z.clear()
            x.append(0)
            y.append(0)
            z.append(0)
        
        #once the photon hits the surface, return the data
        if z[-1] > zmax:
            return x,y,z,theta,phi

#Main function (start point of program --- dummy variables for now)
def main():
    Nonrel_WhiteDwarf(1, 2)
    rel_WhiteDwarf(2, 1)
    plot(1,1, xlabel, ylabel)
    bin_data(0, 0)
    pick_direction()
    pick_optical_depth()       
    move_photon(2, 1)
    
if __name__ == "__main__":
    main()
    

    
