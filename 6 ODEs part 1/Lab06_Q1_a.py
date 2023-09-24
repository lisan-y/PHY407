#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 11:01:51 2021

@author: emmajarvis

Lab06 Q1a
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
    
x_dot = np.linspace(0, 10, 100)

"""
Q1ai
"""
def calc_abs_friction(alpha, beta, vf):
    """
    Calculate the total absolute friction force and plot it as a function of
    velocity.
    
    Parameters
    ----------
    alpha : int
        alpha parameter.
    beta : int
        beta parameter.
    vf : int
        friction speed paramter.

    Returns
    -------
    abs_friction : array
        absolute total friction values to plot.

    """

    Ff = -alpha*x_dot #fluid friction
    Fs = -beta*np.exp(-x_dot/vf) #pseudo static friction

    abs_friction = np.abs(Ff+Fs) #total friction
    
    return abs_friction

#Plotting info
lw = 3
alpha = 0.8

#Plotting for different values of alpha beta and v_f
plt.plot(x_dot, calc_abs_friction(1, 1, 1), color = 'black', 
         label = r'$\alpha = 1, \beta = 1, v_f = 1$', lw = lw, alpha = alpha)

plt.plot(x_dot, calc_abs_friction(2, 5, 2), color = 'tomato', 
         label = r'$\alpha = 2, \beta = 5, v_f = 2$', lw = lw, alpha = alpha)

plt.plot(x_dot, calc_abs_friction(1, 6, 3), color = 'mediumpurple', 
         label = r'$\alpha = 1, \beta = 6, v_f = 3$', lw = lw, alpha = alpha)

plt.plot(x_dot, calc_abs_friction(2, 12, 3), color = 'lightseagreen', 
         label = r'$\alpha = 2, \beta = 12, v_f = 3$', lw = lw, alpha = alpha)

griddy()
plt.xlabel(r'$\dot{x}$ [m/s]')
plt.ylabel('Total Absolute Frictional Force [N]')
plt.legend()
# plt.savefig('/Users/emmajarvis/Desktop/Total_Friction_Force.png', dpi = 300)
plt.show()


