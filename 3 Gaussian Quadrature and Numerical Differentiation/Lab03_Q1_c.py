#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 19:46:54 2021

@author: emmajarvis
"""

import numpy as np
import matplotlib.pyplot as plt

# Plotting parameters
plt.rcParams['figure.figsize'] = (10,6)
plt.rcParams['font.family']='serif'

plt.rcParams['mathtext.fontset'] = 'cm'

S = 20
L = 20
T = 20

plt.rc('font', size = S)
plt.rc('axes', titlesize = T)
plt.rc('axes', labelsize = S)
plt.rc('xtick', labelsize = S)
plt.rc('ytick', labelsize = S)
plt.rc('legend', fontsize = L)
plt.rc('figure', titlesize = S)

def griddy():
    plt.minorticks_on()
    plt.grid(color='black',
             which='major',
             linestyle=":",
             linewidth='0.4',
             )
    plt.grid(color='black',
             which='minor',
             linestyle=":",
             linewidth='0.2',
             )

##############################################################################

#Taken from Lab03_Q1_b
def  gaussxw(N):
    
    #Initial approximation to roots of the Legendre polynomial
    a = np.linspace(3, 4*N-1, N)/(4*N+2)
    x = np.cos(np.pi*a+1/(8*N*N*np.tan(a)))
    
    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = np.ones(N,float)
        p1 = np.copy(x)
        
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta= max(abs(dx))
        
    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)
    return x,w


def gauss(f, x, w, a, b):
    
    xp = 0.5*(b-a)*x + 0.5*(b+a)
    wp = 0.5*(b-a)*w
    s = 0.0
    for k in range(N):
        s += wp[k]*f(xp[k])
    return s

def C(t):
    return np.cos(1/2*np.pi*t**2)

def S(t):
    return np.sin(1/2*np.pi*t**2)

def I(C, S):
    return 1/8*((2*C+1)**2+(2*S+1)**1)
##############################################################################



N=50  # number of slices

x = np.linspace(-3, 10, 100) #array of x values
z = np.linspace(1, 5, 100) # array of z values

def contour_plot(lam):
    """
    Generate contour plot of intensity for (x,z) points

    Parameters
    ----------
    lam : int
        Wavelength in metres.

    Returns
    -------
    None.

    """
    X, Z = np.meshgrid(x, z) #make meshgrid of x and z values for contour plot
    b = X*np.sqrt(2/(lam*Z))  # Lower bound of integral
    a = np.zeros([100, 100])  # upper bound of integral
    

    #initialize list of all fresnel integral values
    Cs = []
    Ss = []

    #get fresnel integral values for each (x,z) value
    y, w = gaussxw(N)
   
    for i in range(len(x)):
        
        Cs.append([])
        Ss.append([])
        for j in range(len(z)):
            #compute using gaussian quadrature
            Cs[i].append(gauss(C, y, w, a[i][j], b[i][j]))
            Ss[i].append(gauss(S, y, w, a[i][j], b[i][j]))
     
    #turn lists into arrays
    Cs = np.asarray(Cs)
    Ss = np.asarray(Ss)

    #compute intensity values
    I_2D = I(Cs, Ss)

    X, Z = np.meshgrid(x, z)
    #make contour plot
    plt.figure()
    plt.contourf(X, Z, I_2D, levels=30, vmin=-0.2, vmax=4.0, cmap = 'magma')
    plt.colorbar(label = 'Intensity')
    plt.xlabel('x [m]')
    plt.ylabel('z [m]')
    plt.title('Intensity of Diffracted Sound Waves \n with a Wavelength of '+str(lam)+' m')
    # plt.savefig('/Users/emmajarvis/Desktop/intensity_2D_lam_'+str(lam)+'_title.png', dpi = 500)
    plt.show()
    
#contour plot for a wavelength of 1m    
contour_plot(1) 

#contour plot for a wavelength ogf 2m
contour_plot(2)

#contour plot for a wavelength of 4m
contour_plot(4)