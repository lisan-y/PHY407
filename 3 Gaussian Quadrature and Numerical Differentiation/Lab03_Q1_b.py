#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 08:59:27 2021

@author: emmajarvis
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sc

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


def gauss(f, N, a, b):
    """
    Evaluate the  integral f from a to b using Gaussian quadrature with 
    N sample points.

    Parameters
    ----------
    N : int
        Number of sample.
    a : float
        Lower integral bound.
    b : float
        Upper integral bound.

    Returns
    -------
    s : float
        Value of integral.

    """
    x,w = gaussxw(N) # get weight values
    
    xp = 0.5*(b-a)*x + 0.5*(b+a)
    wp = 0.5*(b-a)*w
    s = 0.0
    for k in range(N):
        s += wp[k]*f(xp[k]) #value of integral
    return s

#Q1b i

# Define boundaries and slices
N=50  # number of slices

lam = 1 #m, wavelength
z = 3 #m
x = np.linspace(-5,5,100) #m


b = x*np.sqrt(2/(lam*z))  # Lower bound of integral
a = np.zeros(100)  # upper bound of integral


def C(t):
    """
    Fresnel integral C

    """
    return np.cos(1/2*np.pi*t**2)

def S(t):
    """
    Fresnel integral C

    """
    return np.sin(1/2*np.pi*t**2)

def I(C, S):
    """
    Intensity (I/I0) values

    Parameters
    ----------
    C : array of floats
        Fresnel integral.
    S : array of floats
        Fresnel integral.

    Returns
    -------
    array of floats
        Intensity values.

    """

    return 1/8*((2*C+1)**2+(2*S+1)**1)


#initialize list of all fresnel integral values
Cs = []
Ss = []

#get fresnel integral values for each x
for i in range(len(x)):
    #compute fresnel integrals using gaussian quadrature
    Cs.append(gauss(C, N, a[i], b[i])) 
    Ss.append(gauss(S, N, a[i], b[i]))

#turn lists into arrays   
Cs = np.asarray(Cs)
Ss = np.asarray(Ss)


#compute intensity values
I_G = I(Cs, Ss)

#Compute fresnel integrals and intensity values using scipy
S_scipy, C_scipy = sc.fresnel(b)
I_SP = I(C_scipy, S_scipy)


#plot x vs intensity for each method
plt.figure()
plt.plot(x, I_G, label = 'Using Gaussian \n Quadrature', 
          color = 'cornflowerblue', linewidth = 5, alpha = 0.6)
plt.plot(x, I_SP, label = 'Using scipy', 
          color = 'purple', linestyle = '--', linewidth = 2)
griddy()
plt.xlabel('x [m]')
plt.ylabel('I/I0')
plt.title('Intensity of Diffracted Sound Waves')
plt.legend()
plt.tight_layout()
# plt.savefig('/Users/emmajarvis/Desktop/Intensity.png', dpi = 500)
plt.show()

#calculate relative differnce between gaussian quadrature method and scipy fresnel integrals
rel_dif = np.abs(I_SP-I_G)/I_SP

plt.figure()
plt.plot(x, rel_dif, color = 'black', linestyle = '', marker = 'o')
griddy()
plt.xlabel('x [m]')
plt.ylabel('Relative Difference')
plt.title('Relative Difference between the Intensities \n using Scipy and Gaussian Quadrature')
plt.tight_layout()
# plt.savefig('/Users/emmajarvis/Desktop/Intensity_rel_dif.png', dpi = 500)
plt.show()


#Q1b ii


N = np.linspace(3,50,48) #array of N sample points
N = N.astype(int) #make sure each element is an int

#initialize list of all fresnel integral values
Cs = []
Ss = []

#get fresnel integral values for each x and N value
for num in range(len(N)):
    Cs.append([])
    Ss.append([])
    for i in range(len(x)):
        #compute fresnel integrals using gaussian quadrature
        Cs[num].append(gauss(C, N[num], a[i], b[i]))
        Ss[num].append(gauss(S, N[num], a[i], b[i]))


 
Cs = np.asarray(Cs)
Ss = np.asarray(Ss)


#compute maximum relative error for each N value
rel_err_max = []
for i in range(48): 
    I_G = I(Cs[i], Ss[i])
    rel_dif = np.abs(I_SP-I_G)/I_SP
    rel_err_max.append(max(rel_dif[40:])) # restrict domain to where x>0 (approximately)
   
    
   
# Plot maximum relative error as a function of N    
plt.figure()
plt.plot(N, rel_err_max, linestyle = '', marker = 'o', color = 'black')
plt.xlabel('N')
plt.ylabel('Maximum Relative Difference')
plt.title('Maximum Relative Difference as a Function of N points')
griddy()
# plt.savefig('/Users/emmajarvis/Desktop/max_rel_dif.png', dpi = 500)
plt.show()









