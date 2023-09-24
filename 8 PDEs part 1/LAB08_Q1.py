"""
PHY407 Lab 08
Question 1

Animates temperature distribution in a heat conductor
using Gauss-Seidel and overrelaxation.

Author: Lisa Nasu-Yu, Nov 2021
"""

import numpy as np
from pylab import imshow, gray, show
import matplotlib.pyplot as plt

S = 15
L = 15
T = 15

plt.rc('font', size = S)
plt.rc('axes', titlesize = T)
plt.rc('axes', labelsize = S)
plt.rc('xtick', labelsize = S)
plt.rc('ytick', labelsize = S)
plt.rc('legend', fontsize = L)
plt.rc('figure', titlesize = S)

# Constants
top = 200         # Grid squares on a side
side = 80
target = 1e-6  # Target accuracy
w = 0

# Boundary shape
A = [80, 0]
B = [80, 50]
C = [50, 50]
D = [50, 150]
E = [80, 150]
F = [80, 200]
G = [0, 200]
H = [0, 0]

# Create arrays to hold potential values
phi = np.empty([side+1,top+1],float)
# Set boundary conditions
phi[-1, A[1]+1:B[1]+1] = np.linspace(0, 5, 50)
phi[C[0]:-1, B[1]] = np.linspace(7, 5, 30)
phi[C[0], C[1]+1:D[1]+1] = 7
phi[D[0]:-1, D[1]+1] = np.linspace(7, 5, 30)
phi[-1, E[1]+1:] = np.linspace(5, 0, 50)
phi[1:, -1] = np.linspace(5, 0, 80)
phi[0, :] = 10
phi[1:, 0] = np.linspace(5, 0, 80)

old_phi = np.zeros([side+1, top+1], float)

# Main loop
delta = 1.0
# while delta > target:
for t in range(100):
    # Calculate new values of the potential
    for i in range(side):
            for j in range(top):
                if i > 50:
                    # exclude cutout portion of conductor
                    if (j > 51) and (j < 150):
                        pass
                    else:
                        phi[i,j] = (1 + w) * (phi[i+1,j] + phi[i-1,j] \
                                 + phi[i,j+1] + phi[i,j-1])/4 - w*phi[i,j]
                else:
                    phi[i, j] = (1 + w) * (phi[i + 1, j] + phi[i - 1, j] \
                                           + phi[i, j + 1] + phi[i, j - 1]) / 4 - w * phi[i, j]

    # Make a plot
    plt.clf()
    plt.imshow(phi)
    plt.gray()
    plt.colorbar(label='Celsius')
    plt.title("Temperature Distribution in Heat Conductor")
    plt.xlabel('x [cm]')
    plt.ylabel('y [cm]')
    plt.tight_layout()
    plt.draw()
    plt.pause(0.01)

    # Calculate maximum difference from old values
    delta = np.max(np.abs(phi-old_phi))

    # Swap the two arrays around
    old_phi = phi.copy()

print('Animation Complete')

# Make a plot
plt.figure()
plt.imshow(phi)
plt.gray()
plt.colorbar(label='Celsius')
plt.title("Temperature Distribution in Heat Conductor")
plt.xlabel('x [cm]')
plt.ylabel('y [cm]')
plt.tight_layout()
plt.savefig('w0.pdf')

print('Temperature at x=2.5cm, y=1cm: {} Celsius.'.format(phi[-25, 10]))
