"""
PHY407 Final Project

Author: Lisa Nasu-Yu
"""
import numpy as np
import matplotlib.pyplot as plt
from math import exp
from random import random, seed, randint
from time import time

# plot format settings
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


def deriv2(f, i):
    """
    Calculates second derivative
    :param f: [array] series of floats depicting a function
    :param i: [int] index in domain at which to calculate the derivative
    :return: [float]
    """
    d2 = (f[i+1] - f[i]) / dx
    d1 = (f[i] - f[i-1]) / dx
    return (d2 - d1) / dx


def V(i, k):
    """
    Calculates potential
    :param i: [float] position
    :param k: [float] constant in potential
    :return: [float] potential at x
    """
    return 0.5*k*i**2


def hamiltonian(f, i, v, k):
    """
    Calculates hamiltonian of some wavefunction
    :param f: [array] series of floats depicting the wavefunction
    :param i: [int] index in domain at which to calculate
    :param v: numpy function for potential
    :param k: [float] constant in potential
    :return: [float]
    """
    return -1 * (hbar**2/(2*m)) * deriv2(f, i) + v(x[i], k) * f[i]


def normalize(f):
    """
    Normalizes wavefunction
    :param f: [array] series of floats depicting wavefunction
    :return: [array] normalized wavefunction
    """
    return np.sqrt(np.sum([dx*i**2 for i in f]))


def energy(f, v, k):
    """
    Calculates "energy" of trial wavefunction
    :param f: [array] series of floats depicting trial wavefunction
    :param v: numpy function for potential
    :param k: [float] constant in potential
    :return: [float]
    """
    return np.sum([dx*f[i]*hamiltonian(f, i, v, k) for i in range(1, len(x)-2)])


def sim_annealing(v, k):
    """
    Main function for variational Monte Carlo simulation
    :param psi: [array] series of floats depicting initial trial wavefunction
    :param v: numpy function for potential
    :param k: [float] constant in potential
    :return:
    """
    # Main loop
    t = 0
    T = Tmax

    psi = []
    for i in range(len(k)):
        psi.append(0.5 * c0 * np.ones(len(x)) + np.array([np.sum([c_cos[c]*np.cos(c*xi) for c in range(n)]) +
                                                np.sum([c_sin[c]*np.sin(c*xi) for c in range(n)]) for xi in x]))
    psis = []
    times = []
    energys = []
    coeffs = []

    psis.append(psi)
    times.append(t)
    energyn = np.sum([energy(psi[i], v, k[i]) for i in range(len(k))])
    energys.append(energyn)
    c0s = [c0]*len(k)
    c_sins = [c_sin for x_i in range(len(k))]
    c_coss = [c_cos for x_i in range(len(k))]
    # psi = psi

    start_time = time()
    # set seed
    seed(2)
    start_time = time()
    while np.abs(energyn-4.743)/4.743 > 0.05:
        # Cooling
        t += 1
        T = Tmax * exp(-t / tau)

        # Monte Carlo move with gaussian distribution
        oldpsi = np.copy(psi)

        # generate random integer
        for i in range(len(psi)):
            ntest = randint(1, len(x)-2) #xleft, xright
            frac = np.random.normal(0,1)
            dpsi = 0.5*(np.abs(oldpsi[-1][ntest]))
            psi[i][ntest] += frac*dpsi
            psi[i] /= normalize(psi[i]) # normalize

        old_energy = np.copy(energyn)
        energyn = np.sum([energy(psi[i],v, k[i]) for i in range(len(k))])
        delta = old_energy - energyn

        # If the move is rejected, swap them back again
        # if random() > exp(-delta / T):
        if delta < 0:
            psi = oldpsi
            energyn = old_energy
        if t in [1e3, 1e4, 1e5]:
            psis.append(psi)
            times.append(t)
            energys.append(np.sum([energy(psi[i], v, k[i]) for i in range(len(k))]))
        if t > 1e5:
            break
    timer = time() - start_time
    psis.append(psi)
    times.append(t)
    energys.append(np.sum([energy(psi[i], v, k[i]) for i in range(len(k))]))

    return psis, energys, times, timer, coeffs

hbar = 1
dx = 0.2
m=1
x = np.arange(-5,5, dx)
Tmax = 100.0
Tmin = 1e-5
tau = 1e3
bounds = 2
xleft = np.where(np.abs(x+bounds) < 1e-12)[0][0]
xright = np.where(np.abs(x-bounds) < 1e-12)[0][0]
n = 20

c0 = 1
c_sin = np.ones(n)
c_cos = np.ones(n)

psis, energies, ts, timer = sim_annealing(V, [10, 40])
# psis, energies, ts, timer = sim_annealing([trialx], V, [10])

for d in range(np.shape(psis)[1]):
    plt.figure()
    plt.title("{}".format(d))
    for i in range(len(ts)):
        plt.plot(x, psis[i][d], marker='o', label=(ts[i], round(energies[i], 3)))
    plt.legend()
    # plt.savefig('1d montecarlo.pdf')

plt.figure()
plt.title('Energy')
plt.plot(energies)
# plt.savefig('energies.pdf')

plt.figure()
plt.title("Energy {}".format(energies[-1]))
for n, psix in enumerate(psis[-1]):
    plt.plot(x, psix, marker='o', label=n)
plt.legend()
plt.show()
