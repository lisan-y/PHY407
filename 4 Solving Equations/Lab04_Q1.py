#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 10:49:45 2021

@author: emmajarvis
"""

import numpy as np
from numpy.linalg import solve # for LU decomposition
from numpy.random import rand
from time import time
import matplotlib.pyplot as plt

# Plotting parameters
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
    
    
# From Example 6.1

A = np.array([[ 2, 1, 4, 1 ],
[ 3, 4, -1, -1 ] ,
[ 1, -4, 1, 5 ] ,
[ 2, -2, 1, 3 ]], float)
v = np.array( [ -4, 3, 9, 7 ] ,float)


def GaussElim(A_in, v_in):
    """Implement Gaussian Elimination. This should be non-destructive for input
    arrays, so we will copy A and v to
    temporary variables
    IN:
    A_in, the matrix to pivot and triangularize
    v_in, the RHS vector
    OUT:
    x, the vector solution of A_in x = v_in """
    # copy A and v to temporary variables using copy command
    A = np.copy(A_in)
    v = np.copy(v_in)
    N = len(v)

    for m in range(N):
        # Divide by the diagonal element
        div = A[m, m]
        A[m, :] /= div
        v[m] /= div

        # Now subtract from the lower rows
        for i in range(m+1, N):
            mult = A[i, m]
            A[i, :] -= mult*A[m, :]
            v[i] -= mult*v[m]

    # Backsubstitution
    # create an array of the same type as the input array
    x = np.empty(N, dtype=v.dtype)
    for m in range(N-1, -1, -1):
        x[m] = v[m]
        for i in range(m+1, N):
            x[m] -= A[m, i]*x[i]
    return x

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
    x = np.empty(N,float)
    for m in range(N-1,-1,-1):
        x[m] = v[m]
        for i in range(m+1,N):
            x[m] -= A[m, i] *x[i]
    
    return x

print('Q1a')
print('\n For matrix without 0 as first element:')
print('Using gaussian elimination:', GaussElim(A, v))
print('Using partial pivoting:', PartialPivot(A, v))


A = np.array([[ 0, 1, 4, 1 ],
[ 3, 4, -1, -1 ] ,
[ 1, -4, 1, 5 ] ,
[ 2, -2, 1, 3 ]], float)
v = np.array( [ -4, 3, 9, 7 ] ,float)

print('\n For matrix with 0 as first element:')
print('Using partial pivoting:', PartialPivot(A, v))


print('\n Q1b')

Ns = [5, 10, 20, 40, 80, 160, 320, 640]

times_g = []
times_p = []
times_LU = []

errs_g = []
errs_p = []
errs_LU = []

for N  in Ns:
    v = rand(N)
    A = rand(N, N)
    #Gauss Elimination
    start = time()
    x_g = GaussElim(A,v)
    end = time()
    diff = end-start
    times_g.append(diff)
    v_sol = np.dot(A, x_g)
    err = np.mean(abs(v-v_sol))
    errs_g.append(err)
    print(' \n Solution using gaussian elimination for', N, 
          'random values takes a time of ', diff, 'seconds', 
          'and the error is ', err)
    
    #Partial Pivoting
    start = time()
    x_p = PartialPivot(A,v)
    end = time()
    diff = end-start
    times_p.append(diff)
    v_sol = np.dot(A, x_p)
    err = np.mean(abs(v-v_sol))
    errs_p.append(err)
    print(' \n Solution using partial pivoting for', N, 
          'random values takes a time of ', diff, 'seconds', 
          'and the error is ', err)
    
    #LU decomposition
    start = time()
    x_LU = solve(A,v)
    end = time()
    diff = end-start
    times_LU.append(diff)
    v_sol = np.dot(A, x_LU)
    err = np.mean(abs(v-v_sol))
    errs_LU.append(err)
    print(' \n Solution using LU decomposition for', N, 
          'random values takes a time of ', diff, 'seconds', 
          'and the error is ', err)

plt.figure()    
plt.loglog(Ns, errs_g, color = 'lightseagreen', linestyle = '--', marker = 'o', 
             alpha = 0.7, linewidth = 3, label = 'Gaussian Elimination')
plt.loglog(Ns, errs_p, color = 'tomato', linestyle = '--', marker = 'o',
             alpha = 0.7, linewidth = 3, label = 'Partial Pivoting')
plt.loglog(Ns, errs_LU, color = 'mediumpurple', linestyle = '--', marker = 'o',
             alpha = 0.7, linewidth = 3, label = 'LU Decomposition')
plt.ylabel('Errors')
plt.xlabel('Number of values')
griddy()
plt.legend(title = 'Solving Method')
plt.title('Error in Computing x in Ax = v')
plt.savefig('1b_Errors.png', dpi = 400)
plt.show()

plt.figure()
plt.loglog(Ns, times_g, color = 'lightseagreen', linestyle = '--', marker = 'o', 
             alpha = 0.7, linewidth = 3, label = 'Gaussian Elimination')
plt.loglog(Ns, times_p, color = 'tomato', linestyle = '--', marker = 'o',
             alpha = 0.7, linewidth = 3, label = 'Partial Pivoting')
plt.loglog(Ns, times_LU, color = 'mediumpurple', linestyle = '--', marker = 'o',
             alpha = 0.7, linewidth = 3, label = 'LU Decomposition')
plt.ylabel('Times')
plt.xlabel('Number of values')
griddy()
plt.legend(title = 'Solving Method')
plt.title('Time Taken to Compute x in Ax = v')
plt.savefig('1b_Times.png', dpi = 400)
plt.show()

