"""
PHY407 Lab 11
Question 3
Simulated Annealing Optimization
Author: Lisa Nasu-Yu, Dec 2021, adapted from salesman.py, Newman
"""
from math import sqrt, exp
import numpy as np
from numpy import empty
from random import random, randrange, seed
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

N = 25
R = 0.02
Tmax = 10.0
Tmin = 1e-3


# Function to calculate the magnitude of a vector
def mag(x):
    return sqrt(x[0]**2+x[1]**2)


# Function to calculate the total length of the tour
def distance(x):
    s = 0.0
    for i in range(N):
        s += mag(x[i+1]-x[i])
    return s


# Function for simulated annealing
def sim_annealing(nseed, ntau):
    tau = ntau
    # set seed
    ns1 = 1
    seed(ns1)
    # Choose N city locations and calculate the initial distance
    r = empty([N + 1, 2], float)
    for i in range(N):
        r[i, 0] = random()
        r[i, 1] = random()
    r[N] = r[0]
    D = distance(r)
    # Main loop
    t = 0
    T = Tmax

    # set seed
    ns2 = nseed
    seed(ns2)
    while T > Tmin:

        # Cooling
        t += 1
        T = Tmax * exp(-t / tau)

        # Choose two cities to swap and make sure they are distinct
        i, j = randrange(1, N), randrange(1, N)
        while i == j:
            i, j = randrange(1, N), randrange(1, N)

        # Swap them and calculate the change in distance
        oldD = D
        r[i, 0], r[j, 0] = r[j, 0], r[i, 0]
        r[i, 1], r[j, 1] = r[j, 1], r[i, 1]
        D = distance(r)
        deltaD = D - oldD

        # If the move is rejected, swap them back again
        if random() > exp(-deltaD / T):
            r[i, 0], r[j, 0] = r[j, 0], r[i, 0]
            r[i, 1], r[j, 1] = r[j, 1], r[i, 1]
            D = oldD

    return r, D


# Function to plot path
def plot(x, d, save):
    plt.figure()
    plt.title('Path Distance = '+str(d))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.plot(x[:,0], x[:,1], marker='.')
    plt.savefig('lab11_q3_'+save+'.pdf')
    plt.show()
# Set up the graphics
# display(center=[0.5, 0.5])
# for i in range(N):
#     sphere(pos=r[i],radius=R)
# l = curve(pos=r,radius=R/2)

r1, D1 = sim_annealing(1, 1e4)
r2, D2 = sim_annealing(2, 1e4)
r3, D3 = sim_annealing(3, 1e4)

r4, D4 = sim_annealing(1, 1e3)
r5, D5 = sim_annealing(2, 1e5)

print(D1, D2, D3, D4, D5)

rs = [r1, r2, r3, r4, r5]
ds = [D1, D2, D3, D4, D5]

for i in range(len(rs)):
    plot(rs[i], ds[i], str(i))


