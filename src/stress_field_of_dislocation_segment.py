import numpy as np

def stress_field_of_dislocation_segment(mu, nu, r, bx, bz):
    R = np.linalg.norm(r, axis=1)
    sigma_0 = mu/(4*np.pi)/(1-nu)
    sigma=np.zeros((len(r), 6))
    sigma[:,0] = bx*r[:,1]/R**2*(1+(2*(r[:,0]/R)**2))
    sigma[:,1] = bx*r[:,1]/R**2
    sigma[:,2] = 2*bx*nu*r[:,1]/R**2
    sigma[:,3] = (bx*((nu/R)-(r[:,1]**2/R**3))) - (bz*r[:,0]*(1-nu)/R**2)
    sigma[:,4] = -(bx*r[:,0]*r[:,1]/R**3) + (bz*r[:,1]*(1-nu)/R**2)
    sigma[:,5] = -bx*r[:,0]/R**2*(1-(2*(r[:,1]/R)**2))
    
    return sigma_0*sigma
