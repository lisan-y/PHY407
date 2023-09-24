"""
PHY407 Lab 02
Q2 b
Evaluates n-order Bessel functions of the first kind with the Simpson's rule.
Plots the Bessel functions for n = 0, 3, and 5, as well as their difference from the SciPy
module Bessel function, scipy.special.jv.

Author: Lisa Nasu-Yu, 24 Sept 2021
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import jv
from functions_lab02 import *  # make sure in same folder

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


def bessel(n, x):
    """
    Computes an order n Bessel function of the first kind at x,
    using the Simpson's rule with 1000 slices to evaluate the integral.
    :param n: [int] order of Bessel function
    :param x: [float] value at which the Bessel function is evaluated
    :return: [float] value of the order n Bessel function at x
    """
    def f(phi):
        return np.cos(n * phi - x * np.sin(phi))
    return (1 / np.pi) * simpsons(f, 0, np.pi, 1000)


# plot Bessel functions
plt.figure()
plt.title("Bessel Function")
plt.xlabel("x")
plt.ylabel(r'$J_n$')
for n in [0, 3, 5]:
    plt.plot(np.arange(0, 20, 0.1), [bessel(n, i) for i in np.arange(0, 20, 0.1)], label='n='+str(n))
plt.legend()
plt.savefig('Q2b_bessel.pdf')
plt.show()

# plot difference between own Bessel function and scipy.special.jv
plt.figure()
plt.title("Difference Between Simpson's rule Bessel function\nand scipy.special.jv")
plt.xlabel("x")
plt.ylabel('Difference')
for n in [0, 3, 5]:
    plt.plot(np.arange(0, 20, 0.1), jv(n, (np.arange(0, 20, 0.1))) - [bessel(n, i) for i in np.arange(0, 20, 0.1)],
             label='n='+str(n))
plt.legend()
plt.savefig('Q2b_diff.pdf')
plt.show()
