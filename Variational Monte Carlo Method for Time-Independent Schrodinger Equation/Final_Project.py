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


def sim_annealing(psi, v, k):
    """
    Main function for variational Monte Carlo simulation
    :param psi: [array] series of floats depicting initial trial wavefunction
    :param v: numpy function for potential
    :param k: [float] constant in potential
    :return:
    """
    # Main loop
    t = 0  # tracks step number
    # initiate arrays
    psis = []
    times = []
    energys = []

    # input initial values in arrays
    psis.append(psi)
    times.append(t)
    energyn = np.sum([energy(psi[i],v, k[i]) for i in range(len(k))])
    energys.append(energyn)
    psi = psi

    # set seed
    seed(2)
    start_time = time()
    # run loop until within 1 percent of actual value
    while np.abs(energyn-4.743)/4.743 > 0.01:
        t += 1
        oldpsi = np.copy(psi)

        # generate weighted random integer
        psi_max = np.max(psi[0])
        psi_n = 0
        r3 = 1
        while r3 > psi_n / psi_max:
            r1 = random()
            r1 = int(r1*len(x) - (len(x)-1)/2)
            r2 = random()
            r2 = int(r2 * len(x) - (len(x)-1)/2)
            if r2 != 0:
                psi_n = psi[0][r1]*np.cos(np.arctan(r1/r2)) + psi[1][r2]*np.sin(np.arctan(r1/r2))
            else:
                psi_n = psi[0][r1]
            r3 = random()

        # generate random amount to change psi[ntest] by
        frac = np.random.normal(0, 1)
        dpsi = 0.5 * (np.abs(oldpsi[0][r1]))
        psi[0][r1] += frac * dpsi
        psi[0] /= normalize(psi[0])  # normalize

        frac = np.random.normal(0, 1)
        dpsi = 0.5 * (np.abs(oldpsi[1][r2]))
        psi[1][r2] += frac * dpsi
        psi[1] /= normalize(psi[1])  # normalize

        # # generate random integer (not weighted)
        # for i in range(len(psi)):
        #     ntest = randint(1, len(x)-2) #xleft, xright
        #     frac = np.random.normal(0,1)
        #     dpsi = 0.5*(np.abs(oldpsi[-1][ntest]))
        #     psi[i][ntest] += frac*dpsi
        #     psi[i] /= normalize(psi[i]) # normalize

        # compute energy
        old_energy = np.copy(energyn)
        energyn = np.sum([energy(psi[i], v, k[i]) for i in range(len(k))])
        delta = old_energy - energyn

        # If the move is rejected, swap them back again
        if delta < 0:
            psi = oldpsi
            energyn = old_energy

        # store some values in array
        if t in [1e2, 1e3, 1e4, 1e5]:
            psis.append(psi)
            times.append(t)
            energys.append(np.sum([energy(psi[i], v, k[i]) for i in range(len(k))]))

        # prevent overlooping
        if t > 1e6:
            break
        timer = time() - start_time

    # input last sets into array
    psis.append(psi)
    times.append(t)
    energys.append(np.sum([energy(psi[i], v, k[i]) for i in range(len(k))]))
    return psis, energys, times, timer

hbar = 1
m = 1
w = 1
dx = 0.2
x = np.arange(-5,5, dx)
trialx = np.zeros(len(x))
trialy = np.zeros(len(x))
bounds = 2
xleft = np.where(np.abs(x+bounds) < 1e-12)[0][0]
xright = np.where(np.abs(x-bounds) < 1e-12)[0][0]

# set intial trial functions
# trialx[xleft:xright] = 1.0  # constant function
trialx = np.exp(-x**2)
trialx = trialx/ normalize(trialx)  # normalized trial function

# trialy[xleft:xright] = 1.0
trialy = np.exp(-x**2)
trialy = trialy / normalize(trialy)

psis, energies, ts, timer = sim_annealing([trialx, trialy], V, [10,40])

# plot trial wave function progress
for d in range(np.shape(psis)[1]):
    plt.figure()
    plt.title(r'$\psi_{trial}$ for 2-D Harmonic Oscillator')
    plt.ylabel(r'$\psi_{trial}$')
    plt.xlabel('x [m]')
    for i in range(len(ts)):
        plt.plot(x, psis[i][d], label=(ts[i]))
    plt.legend(title='Step')
    # plt.savefig('{}2d montecarlo_exp2.pdf'.format(d))

# plot energies
plt.figure()
plt.title('Energy')
plt.xlabel('Step')
plt.ylabel(r'Energy [J$\cdot \hbar]$')
plt.plot(ts, energies)
# plt.savefig('2energies_exp2.pdf')

# plot final ground state
component = ['x', 'y']
plt.figure()
plt.title("Approximation of Ground State")
plt.ylabel(r'$\psi_{0}$')
plt.xlabel('x [m]')
for n, psix in enumerate(psis[-1]):
    plt.plot(x, psix, marker='o', label=component[n])
plt.legend()
# plt.savefig('2dfinal_exp2.pdf')
plt.show()

# 3d plot
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
X = x
Y = x
X, Y = np.meshgrid(X, Y)
Z = np.zeros((len(x), len(x)))
for i in range(len(x)):
    for j in range(len(x)):
        theta = np.arctan(np.abs(x[j]/x[i]))
        # Z[i,j] = psis[-1][0][i]*np.cos(theta)
        Z[i,j] = psis[-1][0][i]*np.cos(theta) + psis[-1][1][j]*np.sin(theta)

# Plot the surface.
plt.title('Ground State of 2-D Harmonic Oscillator')
surf = ax.plot_surface(X, Y, Z, cmap=plt.cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.zaxis.set_major_locator(plt.LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# set ticks
ax.set_xticks([-5, -2, 0, 2, 5])
ax.set_yticks([-5, -2, 0, 2, 5])
ax.set_zticks([0, 0.5, 1, 1.5])

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel(r'$\psi$')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
# plt.savefig('3d_exp2.pdf')

plt.show()

