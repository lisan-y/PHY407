#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 09:20:02 2021

@author: emmajarvis

Modified code from Lab08 Q2: Exercise 9.5 of Newman (page 431)
to include
"""

import numpy as np
import matplotlib.pyplot as plt

# %matplotlib qt


plt.rcParams['figure.figsize'] = (10,6)
plt.rcParams['font.family']='serif'

plt.rcParams['mathtext.fontset'] = 'cm'

S = 20
L = 15
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
    

def dst(y):
    N = len(y)
    y2 = np.empty(2*N,float)
    y2[0] = y2[N] = 0.0
    y2[1:N] = y[1:]
    y2[:N:-1] = -y[1:]
    a = -np.imag(np.fft.rfft(y2))[:N]
    a[0] = 0.0

    return a

def idst(a):
    N = len(a)
    c = np.empty(N+1,complex)
    c[0] = c[N] = 0.0
    c[1:N] = -1j*a[1:]
    y = np.fft.irfft(c)[:N]
    y[0] = 0.0

    return y

    
#for psi(x, t = 0)    
l = 1 #m, length of string
d = 0.1 #m
C = 1 #m/s
v = 100 #m/s
sigma = 0.3 #m


h = 1e-7 # Time-step
N = 100  # Number of divisions in grid
L = 0.01 # Thickness in meters
a = L/N # grid spacing

# Create empty arrays for phi and psi
phip = np.empty(N+1, float)
psip = np.empty(N+1, float)
phi = np.empty(N+1, float)
psi = np.empty(N+1, float)

#boundary conditions 
phip[0] = 0
phi[N] = 0


c = h * v**2 / a**2 #constant used to calculate phi


# set x to be an array of N values from 0 to l (the length of the string)
x = np.linspace(0, l, N+1) 

t=0
tend = 1e-4
frame_counter = 0



while t<tend:
    
    # Calculate the new values of psi
    
    if t == 0:
        phip[1:N] = np.zeros(N-1) 
        psip = C * x *(l-x) / l**2 * np.exp( ( -(x - d)**2 ) / (2 * sigma**2) )
    
    else:
        phip[1:N] = phi[1:N] + h*psi[1:N]
        psip[1:N] = psi[1:N] + c*(phi[0:N-1]+phi[2:N+1]-2*phi[1 :N])                 
    
    #switch them around    
    phi,phip = phip,phi
    psi,psip = psip,psi
    
    
    # if frame_counter % 20 == 0:
    #     plt.clf()
    #     plt.plot(x, phi, label = str(round(t/1e-6))+'ms', color = 'black')
    #     plt.legend(title = 'Time (s)')
    #     plt.title('Vertical Displacement of Piano String')
    #     plt.ylabel(r'$\Phi(x)$ (m)')
    #     plt.xlabel(r'$x$ (m)')
    #     griddy()
    #     # plt.tight_layout()
    #     plt.xlim(0,1)
    #     plt.ylim(-3e-6, 4e-6)
    #     # if frame_counter % 100 == 0 or frame_counter == 960:
    #     #     plt.savefig('/Users/emmajarvis/Desktop/Q2b_'+str(round(t/1e-6))+'.png', dpi = 200)
    #     plt.draw()
    #     plt.pause(0.01)
    
    t += h #add the timestep
    # frame_counter += 1
    
#New code
L=1  


x = np.linspace(0, l, N) 

ts = np.array([2, 4, 6, 12, 100])*1e-3

for t in ts:
    phi_arr = np.zeros(N)
    phi_arr[0]=0
    phi_arr[N-1]=0
    psi_arr = np.zeros(N)
    psi_arr[0]=0
    psi_arr[N-1]=0
    omegas = []
    for k in range(1,3):
        omega = (k*np.pi/L)
        
        omegas.append(omega)
        phi_arr += np.sin(k*np.pi*x/L)*(np.cos(omega*t))
        psi_arr += np.sin(k*np.pi*x/L)*(np.sin(omega*t))
        phi_arr[0]=0
        phi_arr[N-1]=0
        psi_arr[0]=0
        psi_arr[N-1]=0
        # phi_arr += np.sin(k*np.pi*x/L)
    omegas = np.array(omegas)
    
    # print(phi_arr[3])
    phi_coeffs = dst(phi_arr)
    psi_coeffs = dst(psi_arr)
    # print(phi_coeffs[3])
    
    DerivativeCoeffs = -omega**-2*(-psi_coeffs*np.cos(omega*t)-phi_coeffs*np.sin(omega*t))
  
    # print(DerivativeCoeffs[3])
    d2fphi_dxphi = idst(DerivativeCoeffs)
    # print(d2fphi_dxphi[3])
    
    plt.clf()
    plt.plot(x,d2fphi_dxphi, color = 'black', label = str(t/1e-3))
    plt.ylim(-0.02,0.05)
    plt.legend(title = 'Time', loc = 'upper left')
    griddy()
    plt.ylabel(r'$\phi(x)$ (m)')
    plt.xlabel(r'$x$ (m)')
    plt.title('Vertical Displacement of Piano String')
    # plt.savefig('/Users/emmajarvis/Desktop/Q1_'+str(round(t/1e-3))+'.png', dpi = 200)
    plt.draw()
    plt.pause(0.01)
        
    
    
    