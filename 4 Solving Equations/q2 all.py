"""
PHY407 Lab 04
Question 2

Computes and plots the probability densities of a wave function solving an
asymmetric quantum well.

Author: Lisa Nasu-Yu, Sept 2021
"""

import numpy as np
from numpy.linalg import solve
import matplotlib.pyplot as plt
import matplotlib as mpl


# plot settings
S = 20
L = 15
T = 20

plt.rc('font', size=S)
plt.rc('axes', titlesize=T)
plt.rc('axes', labelsize=S)
plt.rc('xtick', labelsize=S)
plt.rc('ytick', labelsize=S)
plt.rc('legend', fontsize=L)
plt.rc('figure', titlesize=S)

# set plot fonts. runs slower

# mpl.rcParams['legend.frameon'] = False
# mpl.rcParams['figure.autolayout'] = True
#
# plt.rcParams.update({"text.usetex":True, "font.family": 'sans-serif', 'font.sans-serif': ['Helvetica']})
# plt.rcParams.update({"text.usetex":True, "font.family": 'serif', 'font.serif': ['Palatino']})


def Hmatrix(m, n):
    """
    Compute (m,n)th component of the Hamiltonian Operator
    :param m: [int]
    :param n: [int]
    :return: (m,n)th component of the Hamiltonian Operator
    """
    if m == n:
        return 0.5*a + joule2eV * (np.pi**2 * hbar ** 2 * m) / (2 * M * L ** 2)
    elif (m + n) % 2 == 0:
        return 0
    else:
        return (-1) * (8 * a * m * n) / (np.pi**2 * (m**2 - n**2)**2)


def psi(n, x):
    """
    Sums eigenstates to find wave function of the nth energy level
    :param n: [int] energy level
    :param x: [float] x position
    :return: wave function at x
    """
    psi_n = 0
    for m in range(100):
        psi_n += eigs[n][m] * np.sin(np.pi * (m+1) * x / L)
    return psi_n


def gaussxw(N):
    """
    Calculates integration points x and weights w for Gaussian
    quadrature such that sum_i w[i]*f(x[i]) is the Nth-order
    Gaussian approximation to the integral int_{-1}^1 f(x) dx

    Written by Mark Newman <mejn@umich.edu>, June 4, 2011

    :param N: [int] order of Gaussian approximation
    :return: [array] integration points x and weights w
    """

    # Initial approximation to roots of the Legendre polynomial
    a = np.linspace(3, 4*N-1, N)/(4*N+2)
    x = np.cos(np.pi * a + 1 / (8 * N * N * np.tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta > epsilon:
        p0 = np.ones(N, float)
        p1 = np.copy(x)
        for k in range(1, N):
            p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k+1)
        dp = (N+1) * (p0-x*p1) / (1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2 * (N+1) * (N+1) / (N * N * (1-x*x) * dp * dp)

    return x, w


def integrate_psi2(a, b, x, w, n):
    """
    Calculates integral of square of probability density function by gaussian quadrature
    :param a: [float] lower bound of integral
    :param b: [float] upper bound of integral
    :param x: [array] integration points
    :param w: [array] integration weights
    :param n: [int] state number of wave function
    :return: [float] integral of function
    """
    # Map sample points and weights to the integration domain
    xp = 0.5 * (b - a) * x + 0.5 * (b + a)
    wp = 0.5 * (b - a) * w

    # Integrate to calculate period
    s = 0.0
    for i in range(len(xp)):
        s += wp[i] * (psi(n, xp[i])**2)
    return s


# Q 6.9 b

# define constants
L = 5e-10   # width [m]
a = 10  # [eV]
M = 9.1094e-31  # mass [kg]
hbar = 1.054571817e-34  # [J*s]
joule2eV = 6.241509e18  # joules to eV conversion factor


# Q 6.9 c

# set size of H array
mmax = 10
nmax = 10
H = np.zeros((10, 10))

# compute elements of H
for m in range(1, mmax+1):
    for n in range(1, nmax+1):
        H[m-1, n-1] = Hmatrix(m, n)

eigv10 = np.linalg.eigvalsh(H)
print(eigv10)
print('Ground Energy: {} eV'.format(eigv10[0]))


# Q 6.9 d

# set size of H array
mmax = 100
nmax = 100
H = np.zeros((100, 100))

# compute elements of H
for m in range(1, mmax+1):
    for n in range(1, nmax+1):
        H[m-1, n-1] = Hmatrix(m, n)

eigv, eigs= np.linalg.eigh(H)
print(eigv[:10])
print('Ground Energy: {} eV'.format(eigv[0]))

np.savetxt("eigv.csv", np.array([range(1, 11), eigv10, eigv[:10]]).T, delimiter=' & ', newline=' \\\\\n')

# Q 6.9 e
# calculate normalization factor
x, w = gaussxw(20)
A = integrate_psi2(0, L, x, w, 0)

# plot probability densities of ground and 1st 2 excited states
plt.figure()
plt.title(r'Probability Density')
plt.ylabel(r'$|\psi (x)|^2$')
plt.xlabel('x (m)')
plt.plot(np.linspace(0, L, 100), [(psi(0, x))**2 / A for x in np.linspace(0, L, 100)], label='n=0')
plt.plot(np.linspace(0, L, 100), [(psi(1, x))**2 / A for x in np.linspace(0, L, 100)], label='n=1')
plt.plot(np.linspace(0, L, 100), [(psi(2, x))**2 / A for x in np.linspace(0, L, 100)], label='n=2')
plt.legend()
plt.tight_layout()
# plt.savefig('probability_density.pdf')
plt.show()

