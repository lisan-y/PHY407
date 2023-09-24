#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 13:13:57 2021

@author: emmajarvis
"""

import numpy as np
import matplotlib.pyplot as plt


#Plotting Parameters

plt.rcParams['figure.figsize'] = (12,8)
plt.rcParams['font.family']='serif'
plt.rcParams['mathtext.fontset'] = 'cm'
S = 25
L = 15
T = 25
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
#Taken from our solution of Lab 02 Q2a.


def trapezoid(f, a, b, N):
    """
    Evalutes integrals using trapezoid rule with N slices
    :param f: [NumPy function] integrand
    :param a: [float] lower bound of integration
    :param b: [float] upper bound of integration
    :param N: [int] number of slices
    :return: [float] solution to integral with trapezoid rule for N slices
    """
    h = (b - a) / N  # width of slices
    s = 0.5 * f(a) + 0.5 * f(b)  # end slices
    for k in range(1, N):  # add each interior slice
        s += f(a + k * h)
    return h * s


def simpsons(f, a, b, N):
    """
    Evalutes integrals using Simpson's rule with N slices
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

def f(x):
    
    
    """
    Evaluate the function 4 / (1 + x**2) at x
    :param x: [float] value at which to evaluate function
    :return: [float] value of function at x
    """
    return 4 / (1 + x ** 2)


def df(x):
    """
    Evaluate the derivative of the function 4 / (1 + x**2) at x
    :param x: [float] value at which to evaluate the derivative
    :return: [float] value of the derivative x
    """
    return -8 * x / (1 + x ** 2) ** 2

##############################################################################
# Q1a i


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

def gauss(N,a,b):
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

# Define boundaries and slices
Ns=[8, 16, 32, 64, 128, 256, 512, 1024, 2048]  # list of all numbers of slices
a = 0.0  # Lower bound of integral
b = 1.0  # upper bound of integral

#Obtain value of integral for all N values using Gaussian quadrature
gq_vals = []
for N in Ns:
    s = gauss(N, a, b)
    print("The value with Gaussian quadrature method and {} slices is {}".format(N, s))
    gq_vals.append(s)

#Obtain value of integral for all N values using Trapezoid method
trap_vals = []
for N in Ns:
    s = trapezoid(f, a, b, N)
    print("The value with trapezoidal rule and {} slices is {}".format(N, s))
    trap_vals.append(s)

#Obtain value of integral for all N values using Simpson's Method       
simp_vals = []
for N in Ns:
    s = simpsons(f, a, b, N)
    print("The value with Simpson's rule and {} slices is {}".format(N, s))
    simp_vals.append(s)
    


# Q1a ii

#Obtain relative error values for all N values using Gaussian quadrature
errs_gq = []
for i in range(len(gq_vals)-1):
    errs_gq.append(abs(gq_vals[i+1]-gq_vals[i])) #Subtract I_2N - I_N for all N in Ns
 
    
true_val = np.pi #Exact value of the integral
 
#Obtain relative error values for all N values using Trapezoid Method
errs_trap = []
for i in range(len(gq_vals)):
    errs_trap.append(abs(true_val-trap_vals[i])/true_val) #Relative error calculation


#Obtain relative error values for all N values using Simpson's Rule
errs_simp = []
for i in range(len(gq_vals)):
    errs_simp.append(abs(true_val-simp_vals[i])/true_val) #Relative error calculation



# Make the plot of relative error vs N
plt.figure()
plt.loglog(Ns[:-1], errs_gq, linestyle = '', marker = '^', color = 'black', label = 'Gaussian Quadrature')
plt.loglog(Ns, errs_trap, linestyle = '', marker = '<', color = 'red', label = 'Trapezoid')
plt.loglog(Ns, errs_simp, linestyle = '', marker = '>', color = 'blue', label = "Simpson's")
griddy()
plt.xlabel('N')
plt.ylabel('Relative Error')
plt.title('Relative Error as a function of N')
plt.legend(title = 'Integration Method', loc = 'upper right')
# plt.savefig('/Users/emmajarvis/Desktop/Relative_error_1a.png', dpi = 500)
plt.show()