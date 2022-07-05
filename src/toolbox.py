################################
### Some tools
##################################


import numpy as np

def stokeeinstein(radius,kb=1.3806488e-23 , Tcelsius= 35, viscosity=0.693e-3 ):
    """
    Compute the diffusion coefficient from the radius of the particle using the stokes-einstein relationship.
    """
    
    radius=np.array(radius)*1e-9 # nm to m
    
    Tkelvin=273.25+Tcelsius
    
    D=kb*Tkelvin/(6*np.pi*viscosity*radius)*1e4 # in cm2/s
    
    return list(D)

# functions to estimate the max velocity in the PVS :
def A(t, a, f, phi=0, Rv0=8e-4, h0=2e-4) :
    w=2*np.pi*f
    Av0=np.pi*Rv0**2
    Aast0=np.pi*(Rv0+h0)**2
    A0=Aast0-Av0
    return A0*(1+a*np.sin(w*t+phi*2*np.pi))
    
def dAdt(t, a, f, phi=0, Rv0=8e-4, h0=2e-4) :
    w=2*np.pi*f
    Av0=np.pi*Rv0**2
    Aast0=np.pi*(Rv0+h0)**2
    A0=Aast0-Av0
    return A0*(w*a*np.cos(w*t+phi*2*np.pi))
        
def Q (s, t, a, f, l, phi=0, Rv0=8e-4, h0=2e-4) :
    return -dAdt(t, a, f, phi, Rv0,h0)*(s-l)

def U (s, t, a, f, l,  phi=0, Rv0=8e-4, h0=2e-4) :
    return Q(s, t, a, f, l, phi, Rv0, h0)/A ( t, a, f, phi, Rv0, h0)