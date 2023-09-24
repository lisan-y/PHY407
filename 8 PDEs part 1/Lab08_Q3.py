#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 12:00:48 2021

@author: emmajarvis
"""

import numpy as np
import matplotlib.pyplot as plt

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

#time step and spatial step, given in question
dx = 0.02
dt = 0.005


Lx = 2*np.pi #length, given in question
Nx = Lx/dx
Nx = int(Nx)

Tf = 2 #end time, given in question
Nt = Tf/dt
Nt = int(Nt)

epsilon = 1
beta = epsilon*dt/dx

#for plotting the times
error = dt/100

#times to plot
t1 = 0.0
t2 = 0.5
t3 = 1.0
t4 = 1.5

#initial u and u prime arrays with lenth Nx
up = np.empty(Nx+1, float)
u = np.empty(Nx+1, float)

#grid to plot on
xvals = np.linspace(0, Lx, Nx+1)

#boundary conditions
u[0]=0
u[Nx]=0

#initial time
t=0
x=0
  

#iterate until final time
while t < Tf:
    
    #initialize with Euler step
    up[1:Nx] = u[1:Nx] + beta/2*(u[0:Nx-1]+u[2:Nx+1]-2*u[1 :Nx])
    
    u,up = up,u
    
    if t == 0: #initial condition
        up = np.sin(xvals)
    
    else: #use leapfrog method
        up[1:Nx] = u[1:Nx] - beta/2*( ( u[2:Nx+1] )**2 - ( u[0:Nx-1] )**2 )                
    
    #switch them around    
    u,up = up,u
    
    
         
    #plot times 0, 0.5, 1, 1.5    
    if abs(t - t1) < error:
        
        plt.figure()
        plt.plot(xvals, u, color = 'black', label = str(round(t,1)))
        plt.legend(title = 'Time')
        plt.ylabel(r'$u$')
        plt.xlabel(r'$x$')
        griddy()
        # plt.savefig('/Users/emmajarvis/Desktop/Q3_t=0.png', dpi = 200)
        plt.show()
            
            
    if abs(t - t2) < error:
        plt.figure()
        plt.plot(xvals, u, color = 'black', label = str(round(t,1)))
        plt.legend(title = 'Time')
        plt.ylabel(r'$u$')
        plt.xlabel(r'$x$')
        griddy()
        # plt.savefig('/Users/emmajarvis/Desktop/Q3_t=0.5.png', dpi = 200)
        plt.show()
            
            
    if abs(t - t3) < error:
        plt.figure()
        plt.plot(xvals, u, color = 'black', label = str(round(t,1)))
        plt.legend(title = 'Time')
        plt.ylabel(r'$u$')
        plt.xlabel(r'$x$')
        griddy()
        # plt.savefig('/Users/emmajarvis/Desktop/Q3_t=1.png', dpi = 200)
        plt.show()
        
        
    if abs(t - t4) < error:
        plt.figure()
        plt.plot(xvals, u, color = 'black', label = str(round(t,1)))
        plt.ylabel(r'$u$')
        plt.xlabel(r' $x$')
        plt.legend(title = 'Time')
        griddy()
        # plt.savefig('/Users/emmajarvis/Desktop/Q3_t=1.5.png', dpi = 200)
        plt.show()
    
    t+=dt #increase time 
    
    
    

        

      
        
    