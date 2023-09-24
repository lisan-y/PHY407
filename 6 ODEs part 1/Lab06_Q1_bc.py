#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 18:10:55 2021

@author: emmajarvis

Lab 06 Q1b and c
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

"""
Q1bi and Q1bii
""" 
omega = 1
tau = 1
gamma = 1
v_p = 0
v_f = 1


def f(r, t):
    """
    ODE function, no friction

    Parameters
    ----------
    r : x and y
        array.
    t : time
        array.

    Returns
    -------
    array
        solution to ODE.

    """
    x = r[0]
    y = r[1]
    
    spring = -omega**2*(x-v_p*t) #spring component of ODE
    
    # friction = -1/tau*y-gamma*np.exp(-y/v_f)
    
    fx = y
    fy = spring
    
    
    return np.array([fx, fy], float)



a = 0.0
b = 50
h = 1e0

tpoints = np.arange(a,b,h)
xpoints = []
ypoints = []


r = np.array([0.0,1.0],float)
for t in tpoints:
    xpoints.append(r[0])
    ypoints.append(r[1])
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6
    

#calculate energy (potential + kinetic)
energy1 = 1/2*np.array(xpoints)**2 + 1/2*np.array(ypoints)**2
plt.plot(tpoints, energy1, color = 'tomato', label = 'dt = 1', lw = 3, alpha = 0.7)


#Change time step (h)
a = 0.0
b = 50
h = 0.1
tpoints = np.arange(a,b,h)
xpoints = []
ypoints = []

#rk4 method
r = np.array([0.0,1.0],float)
for t in tpoints:
    xpoints.append(r[0])
    ypoints.append(r[1])
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6

energy2 = 1/2*np.array(xpoints)**2 + 1/2*np.array(ypoints)**2

plt.plot(tpoints, energy2, color = 'mediumpurple', label = 'dt = 0.1', lw = 3, alpha = 0.7)
plt.legend()
plt.title(r'Energy Decay of Spring over Time using'+ '\n'+r' $E = 1/2x^2 + 1/2 \dot{x}^2$')
plt.ylabel('Energy [J]')
plt.xlabel('Time [s]')
griddy()
# plt.savefig('/Users/emmajarvis/Desktop/Energy_Decay.png', dpi = 300)
plt.show()

plt.plot(tpoints, xpoints, color = 'black', linewidth = '3')
griddy()
plt.xlabel('Time [s]')
plt.ylabel('x [m]')
plt.title(r'Solution with no friction, $\omega = 1 s^-1$, $v_p = 0$ m/s,' +'\n'+ r' $x_0 = 0$ m and $v_0 = 1$ m/s')
# plt.savefig('/Users/emmajarvis/Desktop/SHM_solution.png', dpi = 300)
plt.show()


"""
Q1biii
"""

omega = 1
tau = 1
gamma = 0.5
v_p = 0.1
v_f = 0.1


def f(r, t):
    """
    ODE function, with friction

    Parameters
    ----------
    r : x and y
        array.
    t : time
        array.

    Returns
    -------
    array
        solution to ODE.

    """
    x = r[0]
    y = r[1]   
    spring = -omega**2*(x-v_p*t)
    if y >= 0: 
        friction = -1/tau*y-gamma*np.exp(-y/v_f) 
    elif y < 0: #condition if velocity is negative
        friction = -1/tau*y-gamma*np.exp(y/v_f)
    fx = y
    fy = spring + friction
    return np.array([fx, fy], float)



a = 0
b = 6
h = 1e-1

tpoints = np.arange(a,b,h)
xpoints = []
ypoints = []

c = v_p


x0 = (1/tau*c +gamma*np.exp(-c/v_f))/(-omega**2)

#rk4 method
r = np.array([x0, c],float)
for t in tpoints:
    xpoints.append(r[0])
    ypoints.append(r[1])
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6
    
    
plt.plot(tpoints, ypoints, color = 'black', linewidth = '3')
griddy()
plt.xlabel('Time [s]')
plt.ylabel('v [m/s]')
plt.title(r'Velocity Values')
plt.ylim(-0.1, 0.2)
# plt.savefig('/Users/emmajarvis/Desktop/ODE_solution_Q1biii_velocity.png', dpi = 300)
plt.show()

plt.plot(tpoints, xpoints, color = 'black', linewidth = '3')
griddy()
plt.xlabel('Time [s]')
plt.ylabel('x (t) [m]')
plt.title(r'Position Values')
# plt.savefig('/Users/emmajarvis/Desktop/ODE_solution_Q1biii_position.png', dpi = 300)
plt.show()

"""
Q1c
"""

#New values of parameters
omega = 1
tau = 1
gamma = 0.5
v_f = 0.1

#v_p is different
v_p = np.linspace(0.1*v_f*np.log(gamma*tau/v_f), 1.5*v_f*np.log(gamma*tau/v_f),5)


def f(r, t, v_p):
    """
    ODE function, with friction, different values of v_p

    Parameters
    ----------
    r : x and y
        array.
    t : time
        array.
    v_p : push velocity
        float.

    Returns
    -------
    array
        solution to ODE.

    """
    x = r[0]
    y = r[1]   
    spring = -omega**2*(x-v_p*t)
    if y >= 0:
        friction = -1/tau*y-gamma*np.exp(-y/v_f) 
    elif y < 0:
        friction = -1/tau*y-gamma*np.exp(y/v_f)
    fx = y
    fy = spring + friction
    return np.array([fx, fy], float)

#plot solution for each value f v_p
for v in v_p:
    a = 0
    b = 7
    h = 1e-1
    
    tpoints = np.arange(a,b,h)
    xpoints = []
    ypoints = []
    
    #rk4 method
    r = np.array([0.0,0.0],float)
    for t in tpoints:
        xpoints.append(r[0])
        ypoints.append(r[1])
        k1 = h*f(r,t, v)
        k2 = h*f(r+0.5*k1,t+0.5*h, v)
        k3 = h*f(r+0.5*k2,t+0.5*h, v)
        k4 = h*f(r+k3,t+h, v)
        r += (k1+2*k2+2*k3+k4)/6
      
        
    plt.plot(tpoints, xpoints, linewidth = 1, label = str(v))
    griddy()
    plt.xlabel('Time [s]')
    plt.ylabel('x [m]')
    plt.legend(title = r'$v_p$ [m/s]', loc = 'upper left')
    plt.title(r'Solution with friction and range of $v_p$')
plt.tight_layout()
# plt.savefig('/Users/emmajarvis/Desktop/ODE_solution_Q1c.png', dpi = 300)
plt.show()


