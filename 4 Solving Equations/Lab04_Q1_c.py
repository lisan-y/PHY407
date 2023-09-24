#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 09:00:48 2021

@author: emmajarvis
"""

import numpy as np
from numpy.linalg import solve
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

def PartialPivot(A, v):
    N = len(v)
    # Gaussian elimination
    for m in range(N):
        # Divide by the diagonal element
        div = A[m,m]
        if div == 0:
            for j in range(len(v)-m):
                i = np.argmax(abs(A[j,:]))
            A[m,: ], A[i, :] = np.copy(A[i, :]), np.copy(A[m, :])
            v[m], v[i] = np.copy(v[i]), np.copy(v[m])
            div = A[m,m]
            
        A [m, : ] /= div
        v[m] /= div
        # Now subtract from the lower rows
        for i in range(m+1,N):
            mult = A [i ,m]
            A[i,:] -= mult*A[m,:]
            v[i] -= mult*v[m]
    
    # Backsubstitution
    x = np.empty(N, dtype = complex)
    for m in range(N-1,-1,-1):
        x[m] = v[m]
        for i in range(m+1,N):
            x[m] -= A[m, i] *x[i]
    
    return x


R1 = 1e3 #ohm
R3 = 1e3 #ohm
R5 = 1e3 #ohm
R2 = 2e3 #ohm
R4 = 2e3 #ohm
R6 = 2e3 #ohm
C1 = 1e-6 #F
C2 = 0.5e-6 #F
xp = 3 #V
w = 1000 #rad/s



A = np.array([[1/R1+1/R4+w*C1*1j, -w*C1*1j, 0],
              [-w*C1*1j, 1/R2+1/R5+w*C1*1j+w*C2*1j, -w*C2*1j],
              [0, -w*C2*1j, (1/R3+1/R6+w*C2*1j)]])

v = [xp/R1, xp/R2, xp/R3]

x = PartialPivot(A,v)
print('\n Method 1')
print('\n The solution for x1, x2 and x3 in complex numbers is ', x)
print('\n The magnitudes of x1, x2 and x3 are ', abs(x))
print('\n The phases of the coefficients are', np.angle(x, deg = True))

x1 = abs(x[0])
x2 = abs(x[1])
x3 = abs(x[2])

t = np.linspace(0,0.015,100)
V1 = x1*np.exp(w*t*1j)
V2 = x2*np.exp(w*t*1j)
V3 = x3*np.exp(w*t*1j)


plt.plot(t, np.real(V1), label = 'V1', color = 'royalblue', linewidth = '3', alpha = 0.7)
plt.plot(t, np.real(V2), label = 'V2', color = 'palevioletred', linewidth = '3', alpha = 0.7)
plt.plot(t, np.real(V3), label = 'V3', color = 'forestgreen', linewidth = '3', alpha = 0.7)
griddy()
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Voltage (V)')
plt.title('Real Voltages')
plt.savefig('Q1c_voltages.png', dpi = 400)
plt.show()

R6 = R6/1000 #H

A = np.array([[1/R1+1/R4+w*C1*1j, -w*C1*1j, 0],
              [-w*C1*1j, 1/R2+1/R5+w*C1*1j+w*C2*1j, -w*C2*1j],
              [0, -w*C2*1j, (1/R3+1/R6+w*C2*1j)]])

v = [xp/R1, xp/R2, xp/R3]

x = PartialPivot(A,v)

print('R6 = R6/omega')
print('\n The solution in complex numbers is ', x)
print('\n The magnitudes of x1, x2 and x3 are ', abs(x))
print('\n The phases of the coefficients are', np.angle(x, deg = True))

x1 = abs(x[0])
x2 = abs(x[1])
x3 = abs(x[2])

t = np.linspace(0,0.015,100)
V1 = x1*np.exp(w*t*1j)
V2 = x2*np.exp(w*t*1j)
V3 = x3*np.exp(w*t*1j)

plt.plot(t, np.real(V1), label = 'V1', color = 'royalblue', linewidth = '3', alpha = 0.7)
plt.plot(t, np.real(V2), label = 'V2', color = 'palevioletred', linewidth = '3', alpha = 0.7)
plt.plot(t, np.real(V3), label = 'V3', color = 'forestgreen', linewidth = '3', alpha = 0.7)
griddy()
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Voltage (V)')
plt.title('Real Voltages')
# plt.savefig('/Users/emmajarvis/Desktop/Desktop/UofT/4th_year/PHY407/Lab_4/Q1c_voltages_R6=iR6', dpi = 400)

