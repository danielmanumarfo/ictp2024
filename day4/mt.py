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

    top = np.array(top)
    rho = np.array(rho)
     
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

def plotrhophi2(rho,phi,rho2,phi2,period):
    
    plt.figure(figsize=(6,12))
    ##fig, axs = plt.subplots(2,1,figsize=(16,9))
    appres = plt.subplot(211)
    plt.plot(period,rho,'ro-',label = 'model 01')
    plt.plot(period,rho2,'bo-', label = 'model 02')
    plt.xlabel('period [s]')
    plt.ylabel(r'$\rho_a [\Omega m]$')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim((0.1,10000))
    plt.grid()
    #plt.grid(b=True,which='both')
    plt.title('FWD response: App. resistivity & phase')
    plt.legend()
    #plt.set_figheight(2)
    
    phase = plt.subplot(212)
    plt.plot(period,phi,'ro-')
    plt.plot(period,phi2,'bo-')
    plt.xlabel('period [s]')
    plt.ylabel('phase [°]')
    plt.xscale('log')
    plt.grid()
    plt.ylim((0,90))
    #plt.grid(b=True,which='both')

def plotrhophi(rho,phi,period,title = 'FWD response',myfmt = 'ro-'):
    
    plt.figure(figsize=(6,12))
    appres = plt.subplot(211)
    plt.plot(period,rho,myfmt)
    plt.xlabel('period [s]')
    plt.ylabel(r'$\rho_a [\Omega m]$')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    plt.ylim((0.1,10000))
    #plt.grid(b=True,which='both')
    plt.title(title+': App. resistivity & phase')
    
    
    phase = plt.subplot(212)
    
    plt.plot(period,phi,myfmt)
    plt.xlabel('period [s]')
    plt.ylabel('phase [°]')
    plt.xscale('log')
    plt.grid()
    plt.ylim((0,90))
    #plt.grid(b=True,which='both')

def plotrhophi_eb(rho,phi,rho_err,phi_err,rho2,phi2,period,title = 'FWD response'):
    
    fig = plt.figure(figsize=(6,12))
    appres = plt.subplot(211)
    plt.errorbar(period,rho,yerr = rho_err, color = 'black',fmt = 'o',label = 'observed')
    plt.plot(period,rho2,'--r',label = 'modelled')
    #plt.plot(period,rho2,'ro-', label = 'model 02')
    plt.xlabel('period [s]')
    plt.ylabel(r'$\rho_a [\Omega m]$')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim((0.1,10000))
    fig.set_figheight(2)
    plt.grid()
    #plt.grid(b=True,which='both')
    plt.title(title+': App. resistivity & phase')
    plt.legend()
    
    phase = plt.subplot(212)
    plt.errorbar(period,phi,yerr = phi_err, color = 'black',fmt = 'o',label = 'observed')
    plt.plot(period,phi2,'r--')
    plt.xlabel('period [s]')
    plt.ylabel('phase [°]')
    plt.xscale('log')
    plt.grid()
    plt.ylim((0,90))
    #plt.grid(b=True,which='both')

def plotmodel2(res,thickness,res2,thickness2,tol = '--r', tol2='--b'):
    
    # convert the thicknesses in a depth vector from 0 to after the maximum depth
    depth = np.zeros(len(thickness)+2)
    depth[1:-1] = np.cumsum(thickness[:])
    depth[-1] = np.sum(thickness)+500
    # resampling the resistivity vector so it includes the first layer from depth 0 m
    rho = np.zeros(len(res)+1)
    rho[0:-1] = res[0:]
    rho[-1] = res[-1]
    # convert the thicknesses in a depth vector from 0 to after the maximum depth
    depth2 = np.zeros(len(thickness2)+2)
    depth2[1:-1] = np.cumsum(thickness2[:])
    depth2[-1] = np.sum(thickness2)+500
    # resampling the resistivity vector so it includes the first layer from depth 0 m
    rho2 = np.zeros(len(res2)+1)
    rho2[0:-1] = res2[0:]
    rho2[-1] = res2[-1]
    #starting the figure plot 
    plt.figure()
    plt.step(rho,depth,tol, label = 'model 01')
    plt.step(rho2,depth2,tol2, label = 'model 02')
    plt.title('rho vs depth')
    plt.xlabel('rho [ohm.m]')
    plt.ylabel('depth [m]')
    plt.xscale('log')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.xlim((0.01,10000))
    
    plt.show()

def plotmodel(res,thickness,tol='--b'):
    #tol = '--b'
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
    plt.xlim((0.01,10000))
    
    
    
    plt.show()

def give_solution():
    # model
    rho = [80,10,2.68,10, 300] # resistivity of each layer in [Ohmm]
    # --- DO NOT EDIT BELOW --- #
    thickness =[10, 600, 1000, 500] # thickness of layers in [m], separated by commas, last layer is assumed to have infinite thickness --> thickness has one entry less than rho
    # read data from textfile
    data = np.genfromtxt('data.dat')
    per = 1./data[4,:]
    rhoa, phi = waitMT(thickness,rho,per) # FWD response of model 

    return rhoa, phi
