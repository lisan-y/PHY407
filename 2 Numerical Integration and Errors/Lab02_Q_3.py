#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 11:06:02 2021

@author: emmajarvis
"""

import numpy as np
from scipy.constants import epsilon_0
from scipy.special import kn
import matplotlib.pyplot as plt

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
plt.rcParams['figure.figsize'] = (10,6)
plt.rcParams["font.family"] = "serif"

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

"""
Function simpsons taken from functions_lab02.py
"""
def simpsons(f, a, b, N):
    """
    Evalutes integrals using Simpson's rule with N slices using the formula
    I(a,b) = 1/3 * h * [f(a) + f(b) + 4 * sum(f(a + k*h)) + 2 * sum(f(a + m*h))]
    for odd k in [1, N-1] and even m in [2, N-2],
    derived from manipulation of the Taylor expansion of f.
    :param f: [NumPy function] integrand
    :param a: [float] lower bound of integration
    :param b: [float] upper bound of integration
    :param N: [int] number of slices
    :return: [float] solution to integral with Simpson's rule for N slices
    """
    h = (b - a) / N  # width of slices
    s = f(a) + f(b)  # end slices
    for k in range(1, N, 2):  # odd terms
        s += 4 * f(a + k * h)
    for k in range(2, N, 2):  # even terms
        s += 2 * f(a + k * h)
    return (h / 3) * s


#3a

Q = 1e-13 #C
l = 1/1000 #m
z = 0 #m
r = np.linspace(0.25/1000, 5/1000, 100) #m
k0 = kn(0, r**2/(2*l**2)) #zeroth order Bessel function of the second kind

"""
FUNCTION: V

Purpose: compute the electrostatic potential integral for the simpson' rule
"""
def V(u):
    num = Q*np.exp(-(np.tan(u))**2)
    denom = 4*np.pi*epsilon_0*np.cos(u)**2*np.sqrt((z-l*np.tan(u))**2+r**2)
    return num/denom

"""
FUNCTION: V_known

Purpose: 
    To compute the electrostatic potential of the known solution with z = 0.

Paramter:
    r: array of radii

Output:
    Value of the integral
"""
def V_known(r):
    num = Q*np.exp(r**2/(2*l**2))*k0
    denom = 4*np.pi*epsilon_0*l
    return num/denom

    
integral = simpsons(V, -np.pi/2, np.pi/2, 8)
known = V_known(r)

plt.plot(r, known, label = 'Known Result', 
          color = 'deepskyblue', alpha = 1, linewidth = 3,linestyle = '--')
plt.plot(r, integral, label = "Calculated using Simpson's Rule \n with N = 8 slices", 
          color = 'crimson', alpha = 0.7, linewidth = 3,linestyle = '-')

plt.xlabel('Radius [m]')
plt.ylabel('Potential [V]')
plt.title('Electrostatic Potential for a Line Charge')
plt.legend()
griddy()
plt.tight_layout()
# plt.savefig('potentials.png', dpi=500)
plt.show()


integral = simpsons(V, -np.pi/2, np.pi/2, 50)
dif = abs(known - integral)/known

plt.plot(r, dif, color = 'darkviolet', 
          linewidth = 2, linestyle = '-')
plt.xlabel('Radius [m]')
plt.ylabel('Difference')
plt.title("Difference between the known result \n and using the Simpson's Method")
griddy()
plt.tight_layout()
# plt.savefig('difference.png', dpi=500)
plt.show()



# 3.b
z = np.linspace(-5/1000, 5/1000, 100) #m
r = np.linspace(0.25/1000, 5/1000, 100) #m
X, Y = np.meshgrid(r, z) #prepare for contour and surface plots


"""
New definition of V to use an array of values for z instead of just 0.
"""
def V_3b(u):
    num = Q*np.exp(-(np.tan(u))**2)
    denom = 4*np.pi*epsilon_0*np.cos(u)**2*np.sqrt((Y-l*np.tan(u))**2+X**2)
    return num/denom


# 3D Surface Plot
V = simpsons(V_3b, -np.pi/2, np.pi/2, 50)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(X, Y, V, cmap = 'magma',
                       linewidth = 0, antialiased = True)
fig.colorbar(surf, shrink=0.5, aspect=8, label = 'Potential [V]')
ax.set_xlabel('Radius [m]')
ax.set_ylabel('z-axis [m]')
plt.show()

# Contour Plot
plt.contourf(Y, X, V, 20, cmap='magma')
plt.colorbar(label = 'Potential [V]');
plt.xlabel('z axis [m]')
plt.ylabel('Radius [m]')
plt.title('Potential Field')
plt.tight_layout()

# plt.savefig('Potential_field_3d.png', dpi=500)
 




