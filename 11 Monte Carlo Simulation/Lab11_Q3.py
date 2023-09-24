"""
PHY407 Lab 11
Question 3

Author: Lisa Nasu-Yu, Dec 2021, adapted from salesman.py, Newman
"""
from math import sqrt, exp
import numpy as np
from random import random, randrange, seed
import matplotlib.pyplot as plt
from time import time

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

N = 25
R = 0.02
Tmax = 10.0
Tmin = 1e-3
tau = 1e4


# Computes function defined in Q3bi
def f1(i, j):
    return i**2 - np.cos(4*np.pi*i) + (j-1)**2


# Computes function defined in Q3bii
def f2(i, j):
    return np.cos(i) + np.cos(np.sqrt(2)*i) + np.cos(np.sqrt(3)*i) + (j-1)**2


# Function for simulated annealing
def sim_annealing(function):
    # Main loop
    t = 0
    T = Tmax

    # initial position
    x = 2
    y = 2
    # initiate x and y arrays
    xarr = [x]
    yarr = [y]
    tarr = [0]
    start_time = time()
    # set seed
    ns1 = 2
    seed(ns1)
    while T > Tmin:

        # Cooling
        t += 1
        T = Tmax * exp(-t / tau)

        # Monte Carlo move with gaussian distribution
        oldx, oldy = x, y
        dx, dy = np.random.normal(0, 1, 2)
        x, y = x+dx, y+dy

        delta = function(x, y) - function(oldx, oldy)

        # If the move is rejected, swap them back again
        if random() > exp(-delta / T):
                x = oldx
                y = oldy

        # FOR 3bii
        elif (0 < x < 50) and (-20 < y < 20):
            xarr.append(x)
            yarr.append(y)
            tarr.append(time()-start_time)
        else:
            x = oldx
            y = oldy

        # FOR 3bi
        # else:
        #     xarr.append(x)
        #     yarr.append(y)
        #     tarr.append(time()-start_time)

    return xarr, yarr, tarr


# Function to plot path
def plot(x, y, ylabel, save):
    plt.figure()
    plt.title('Final {} = {}'.format(ylabel, y[-1]))
    plt.xlabel('Time [s]')
    plt.ylabel(ylabel)
    plt.plot(x, y, '.')
    # plt.savefig('lab11_q3_'+save+'.pdf')
    plt.show()


# # Q3bi
# xbi, ybi, tbi = sim_annealing(f1)
# plot(tbi, xbi, 'x', 'bix')
# plot(tbi, ybi, 'y', 'biy')

# Q3bii
xbii, ybii, tbii = sim_annealing(f2)
plot(tbii, xbii, 'x', 'biix')
plot(tbii, ybii, 'y', 'biiy')




