
import numpy as np
import matplotlib.pyplot as plt 


#Random generator from 0 to 1
rng = np.random.default_rng()


def pick_direction():
    theta = np.arccos(1 - 2 * rng.random())
    phi = 2*np.pi*rng.random()
    return theta, phi

def test_dir(n):
    arrtheta = np.zeros(n)
    arrphi = np.zeros(n)
    
    for i in range(n):
        arrtheta[i], arrphi[i] = pick_direction()
    
    x = np.sin(arrtheta)*np.cos(arrphi)
    y = np.sin(arrtheta)*np.sin(arrphi)
    z = np.cos(arrtheta)
    

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x,y,z, marker='.')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()
      
def pick_optical_depth():
    x = rng.random()
    tau = - np.log(1-x)
    return tau

def test_op_depth(n,nbin):
    arrdepth = np.zeros(n)
    
    for i in range(n):
        arrdepth[i] = pick_optical_depth()
    
    dmax = np.amax(arrdepth) 
    dmin = np.amin(arrdepth)
    
    bins = np.zeros(nbin)
    
    for i in range(n):
        b = int((arrdepth[i] - dmin)/(dmax+1e-7 - dmin) * nbin)
        bins[b] += 1
        
    plt.bar(range(nbin), bins)
    plt.show()
    


def move_photon(tmax, zmax):
    x = []
    y = []
    z = []
    
    x.append(0)
    y.append(0)
    z.append(0)
    
    while True:
        theta, phi = pick_direction()
        opdepth = pick_optical_depth()
        
        s = opdepth/tmax
        dx = s * np.sin(theta)*np.cos(phi)
        dy = s * np.sin(theta)*np.sin(phi)
        dz = s * np.cos(theta)
        
        x.append(x[-1] + dx)
        y.append(y[-1] + dy)
        z.append(z[-1] + dz)
        
        print( f"Photon Position: {x[-1]},{y[-1]},{z[-1]}" )
        
        if z[-1] < 0:
            x.clear()
            y.clear()
            z.clear()
            x.append(0)
            y.append(0)
            z.append(0)
        
        if z[-1] > zmax:
            return x,y,z,theta, phi
        
x,y,z,theta,phi = move_photon(10,1)
print(theta, phi)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x,y,z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()

        
        
    
    