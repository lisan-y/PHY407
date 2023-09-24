"""
PHY407 Lab 06
Question 2

Solves the system of an N-storey vibrating building using the Verlet method
and the normal modes.

Author: Lisa Nasu-Yu, Oct 2021
"""

import numpy as np
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


def A_matrix(N):
    """Create matrix for N dimensions, as defined in lab assignment
    :param N: [int] number of storeys
    :return: [array] matrix A
    """
    # first row
    A = np.zeros((N, N))
    A[0, 0] = -2
    A[0, 1] = 1
    # middle rows
    for n, row in enumerate(A[1:-1]):
        A[n+1, n] = 1
        A[n + 1, n + 1] = -2
        A[n + 1, n + 2] = 1
    # last row
    A[-1, -1] = -2
    A[-1, -2] = 1
    return k2m * A


def verlet(N, x0, A, title, eigen=False, valn=0):
    """
    Solves and plots the system of vibrating building by Verlet method
    :param N: [int] number of storeys
    :param x0: [array] initial displacements of each floor
    :param A: [array] matrix describing system
    :param valn: [int] normal mode number (if computing normal modes)
    :param title: [str] title of plot
    :param eigen: [bool] True if simulating normal modes
    :return: time series plot of motion of the system
    """
    A = A
    # initial conditions
    x = x0
    v = np.zeros(N)  # all floors start at rest

    # initiate arrays
    xpoints = np.zeros((len(tpoints), N))
    vpoints = np.zeros((len(tpoints), N))

    # Verlet method computation
    for t in range(len(tpoints)):
        xpoints[t] = x
        vpoints[t] = v
        x += dt * v
        v += dt * np.dot(A, x)

    plt.figure()
    plt.title(title)
    for n in range(N):
        plt.ylabel('Displacement [m]')
        plt.xlabel('Time [s]')
        if eigen:
            # plot normal modes
            plt.plot(tpoints, np.array(xpoints)[:, n] * np.cos(np.sqrt(-eigval[valn]) * tpoints), label=n)
        else:
            plt.plot(tpoints, np.array(xpoints)[:, n], label=n)
    plt.legend(title='Floor', loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    # plt.savefig(title+'.pdf')
    plt.show()


# Q2 a)

# Define constants
k2m = 400  # [rad/s^2]
dt = 1e-3  # [s]
tpoints = np.arange(0, 1, dt)

# solve 3 and 10 storey building
for N in [3, 10]:
    # initial positions
    x0 = np.zeros(N)
    x0[0] = 0.1  # [m] set first floor at 10cm, all others at 0m
    verlet(N, x0, A_matrix(N), str(N)+'-Storey Oscillating Building')


# Q2 b)

N = 3  # number of stories
A = A_matrix(N)

# compute eigenvalues and vectors
eigval, eigvec = np.linalg.eigh(A)

# solve system for normal modes
for n, x0 in enumerate(eigvec):
    verlet(N, x0, eigval, 'Normal Mode '+str(n), eigen=True, valn=n)
