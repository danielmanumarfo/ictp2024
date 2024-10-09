# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 17:44:31 2022

@author: jcruces
"""

# ==============================================================================
import numpy as np
import matplotlib.pyplot as plt
# import os
# ==============================================================================



def waitMT(top,rho,p):
     
    dl = top
    nper = p.size # no of periods
    nl = top.size # no of layers

    if (top.size+1 != rho.size):
        print("Number of layer thicknesses and resistivities inconsistent. No FWD calculation.\n")
        return
    else:
        print(f"Calculating 1D FWD ...")

        
    

    z = np.empty(nper,dtype=complex)
    
     
    u=4*10**-7*np.pi  #\mu_0
    
    w = 2*np.pi/p #use frequencies
    
    ## calculate c response recursively:
    
    for ip in range(0,nper): # loop over periods
        g=complex(1,0) 
        
        for il in range (nl-1,-1,-1): #loop over layers, starting from bottom
            a = np.tanh(np.sqrt(complex(0,u*w[ip]) / rho[il])*dl[il])
            
            k=g*np.sqrt(rho[il+1]/rho[il])
            
            g=(k+a)/(1+k*a)
       
        #Calculate and store the impedances
        z[ip] = complex(0,w[ip])*0.001*g/np.sqrt(complex(0,u*w[ip]/rho[0])); # z response
        
    rhoa, phi = rhophi(z,p)
    
    return rhoa, phi

def rhophi(z,p):
     
    #deg = 180/np.pi
    rhoa = 0.2*p*np.abs(z)**2
    phi = np.angle(z,deg=True)
    
    return rhoa, phi

def plotrhophi(rho,phi,period):
    
    plt.figure()
    appres = plt.subplot(211)
    plt.plot(period,rho,'ro-')
    plt.xlabel('log10(period) [s]')
    plt.ylabel(r'$\rho_a [\Omega m]$')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(b=True,which='both')
    
    phase = plt.subplot(212)
    plt.plot(period,phi,'ro-')
    plt.xlabel('log10(period) [s]')
    plt.ylabel('phase [°]')
    plt.xscale('log')
    plt.grid()
    plt.grid(b=True,which='both')


def plotmodel(res,thickness):
    tol = '--b'
    # convert the thicknesses in a depth vector from 0 to after the maximum depth
    depth = np.zeros(len(thickness)+2)
    depth[1:-1] = np.cumsum(thickness[:])
    depth[-1] = np.sum(thickness)+500
    # resampling the resistivity vector so it includes the first layer from depth 0 m
    rho = np.zeros(len(res)+1)
    rho[0:-1] = res[0:]
    rho[-1] = res[-1]
    #starting the figure plot 
    plt.figure()
    plt.step(rho,depth,tol)
    plt.title('rho vs depth')
    plt.xlabel('rho [ohm.m]')
    plt.ylabel('depth [m]')
    plt.xscale('log')
    plt.gca().invert_yaxis()
    
    
    plt.show()