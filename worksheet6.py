import numpy as np
import matplotlib.pyplot as plt

def solveODE(x_start, x_stop, h, y0, z0, f, g, stop = lambda x, y, z, data: False, data = None):
    """
    	def solveODE(x_start, x_stop, h, y0, z0, f, g):
    	
        Solves an ODE using RK4
        
        Parameters: 
          * x_start - starting point
          * x_stop - stopping
          * h - step size
          * y0, z0 - initial condition
          * f, g - ODE functions
          
        Returns:
          * x, y, z - arrays that hold the solution data
    """
   
    
    #initialize the arrays for each variable
    x = np.arange(x_start,x_stop,h)
    N = len(x)
    y = np.zeros(N)
    z = np.zeros(N)
    
    #set initial conditions
    y[0] = y0
    z[0] = z0
    
    # default to runge-kutta
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
        
        if stop(x[i+1], y[i+1], z[i+1], data):
        	return x[0:i], y[0:i], z[0:i]
        
    # return arrays of x,y,z variables
    return x,y,z   

def Nonrel_WhiteDwarf(x_stop, rho_c = 4.045e6, error=1e-4):
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
        
    def stop(eta, rho, sigma, data):
        if abs(rho)  < data:
            print(f"### Stopping the ODE solver at radius eta = {eta}")
            return True
        else:
            return False
        
    #call solveODE to solve for values of eta, rho, sigma
    eta,rho,sigma = solveODE(0.000001,x_stop,0.001,1,0,f,g, stop, error)

       
    
   
    
    #write formula for lambda then convert to km
    lamb = (((n + 1) * k * rho_c**((1-n)/n))/ (4 * np.pi * G))**(0.5) #in cm
    lamb = lamb * 10**(-5)  #convert to km
    
    print(f"Lambda = {lamb:.1f} km")
    
    #calculate dimensionless mass
    y = rho**n * eta**2
    m = np.trapz(y,eta)
    
    print(f"Dimensionless mass is {m:.2f}")
    
    #Calculate mass M
    #within the formula multiply lambda by 10^5 to convert to cm
    # our final asnwer for M will be in grams, so we will multiply the entire 
    #formula by 10^-3 to convert to kg
    
    M = 4 * np.pi * rho_c * (lamb * 10**(5))**3 * m * 10**(-3)
    print(f"Physical mass is {M:.2e} kg")
    
    #find real density and radius of the star
    eta = eta * lamb  # convert to km using the constant lambda
    rho = rho * rho_c 
    radius = eta[-1]
     
    
        
    return eta, rho, M, radius

def main():

    eta, rho, M, r = Nonrel_WhiteDwarf(4.0, rho_c = 4.045e7)
    plt.plot(eta, rho)
    plt.show()
    
if __name__ == "__main__":
    main()
