"""
PHY407 Lab 05
Question 2 a) Re-visiting the relativistic spring

Simulates the motion of a relativistic spring using the Euler-Cromer method
and by Fourier Transform.

Computes the estimated period

Author: Lisa Nasu-Yu, Oct 2021
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


# plot settings
S = 15
L = 15
T = 15

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


# Functions from Lab 3
############################


def g(x, x0):
    """
    Computes the speed of a particle at position x, with initial position x0, mass m, and spring constant k
    :param x: [float] position at which g(x) is evaluated
    :param x0: [float] initial position of particle
    :return: [float] speed of particle at x
    """
    return c * np.sqrt((k * (x0**2 - x**2) * (2*m*c**2 + k*(x0**2 - x**2)/2)) /
                       (2 * (m*c**2 + k*(x0**2 - x**2)/2) ** 2))


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


def period(x0, x, w):

    # Define integral bounds
    a = 0
    b = x0

    # Map sample points and weights to the integration domain
    xp = 0.5 * (b - a) * x + 0.5 * (b + a)
    wp = 0.5 * (b - a) * w

    # Integrate to calculate period
    s = 0.0
    for i in range(len(x)):
        s += wp[i] * (g(xp[i], x0) ** -1)
    return 4 * s, xp, wp

##################################################


def euler_cromer(x0, t_f, dt):
    """
    Simulates the position and velocity of a relativistic spring by Euler-Cromer method.
    :param x0: [float] initial position, metres
    :param t_f: [float] length of simulation, seconds
    :param dt: [float] time step
    :return: positions and velocities through the time t_f
    """
    # Initiate time, position, and velocity arrays
    t = np.arange(0, t_f, dt)  # time [s]
    x = np.zeros(len(t))  # position [m]
    v = np.zeros(len(t))  # speed [m/s]

    # Set initial conditions in arrays
    x[0] = x0  # [m]
    v[0] = 0  # [m/s]

    # Calculate position & velocity in increments
    for i in range(len(t) - 1):
        acceleration = -1 * (k / m) * x[i] * (1 - (v[i]**2 / c**2)) ** (3/2)  # [m/s**2]
        v[i + 1] = v[i] + acceleration * dt
        x[i + 1] = x[i] + v[i + 1] * dt

    return x, v


def plotxv(x, v, t, save):
    """
    Plots position and speed of a relativistic spring
    :param x: [array] position
    :param v: [array] speed
    :param t: [float] total time, seconds
    :param save: [str] name of file to save
    :return: plot of position and speed
    """
    plt.figure()
    plt.title(r'Relativistic Spring, $x_0={}$'.format(save))
    plt.xlabel('Time [s]')
    plt.plot(np.arange(0, t, dt), x, label='Position [m]')
    plt.plot(np.arange(0, t, dt), v, label='Speed [m/s]')
    plt.legend()
    plt.savefig('spring'+save+'.pdf')
    plt.show()


# Define constants
c = 3e8  # speed of light [m/s]
k = 12  # spring constant [N/m]
m = 1  # mass [kg]

# Define initial positions [m]
x01 = 1
x02 = c / (2 * np.sqrt(3))  # x_c, initial position to reach v=c, calculated in lab 3
x03 = x02 * 10

# Define length of each simulation
t1 = 20
t2 = 25
t3 = 150

# Define time step
dt = 0.0001

x1, v1 = euler_cromer(x01, t1, dt)
x2, v2 = euler_cromer(x02, t2, dt)
x3, v3 = euler_cromer(x03, t3, dt)

# plot positions and velocities
plotxv(x1, v1, t1, '1')
plotxv(x2, v2, t2, 'x_c')
plotxv(x3, v3, t3, '10x_c')

# calculate period with gaussian quadrature, N=200
x, w = gaussxw(400)
T1, _, _ = period(x01, x, w)
T2, _, _ = period(x02, x, w)
T3, _, _ = period(x03, x, w)

x1ft = np.fft.rfft(x1)
x2ft = np.fft.rfft(x2)
x3ft = np.fft.rfft(x3)

v1ft = np.fft.rfft(v1)
v2ft = np.fft.rfft(v2)
v3ft = np.fft.rfft(v3)

# Compute frequencies of fft
freq1 = np.fft.rfftfreq(len(x1), dt)
freq2 = np.fft.rfftfreq(len(x2), dt)
freq3 = np.fft.rfftfreq(len(x3), dt)

plt.figure()
plt.title('Normalized Position')
plt.xscale('log')
plt.xlabel('w [rad/s]')
plt.ylabel(r'$|\hat{x}(w)|/|\hat{x}(w)|_{max}$')
plt.vlines(2 * np.pi/T1, 0, 1, color='blue')
plt.vlines(2 * np.pi/T2, 0, 1, color='orange')
plt.vlines(2 * np.pi/T3, 0, 1, color='green')
plt.plot(freq1 * 2 * np.pi, np.abs(x1ft) / np.max(np.abs(x1ft)), color='blue', label='1')
plt.plot(freq2 * 2 * np.pi, np.abs(x2ft) / np.max(np.abs(x2ft)), color='orange', label=r'$x_c$')
plt.plot(freq3 * 2 * np.pi, np.abs(x3ft) / np.max(np.abs(x3ft)), color='green', label=r'$10x_c$')
plt.legend(title='Initial Position [m]')
plt.tight_layout()
plt.savefig('springftx.pdf')
plt.show()

plt.figure()
plt.title('Normalized Velocity')
plt.xscale('log')
plt.xlabel('w [rad/s]')
plt.ylabel(r'$|\hat{v}(w)|/|\hat{v}(w)|_{max}$')
plt.vlines(2 * np.pi/T1, 0, 1, color='blue')
plt.vlines(2 * np.pi/T2, 0, 1, color='orange')
plt.vlines(2 * np.pi/T3, 0, 1, color='green')
plt.plot(freq1 * 2 * np.pi, np.abs(v1ft) / np.max(np.abs(v1ft)), color='blue', label=r'1')
plt.plot(freq2 * 2 * np.pi, np.abs(v2ft) / np.max(np.abs(v2ft)), color='orange', label=r'$x_c$')
plt.plot(freq3 * 2 * np.pi, np.abs(v3ft) / np.max(np.abs(v3ft)), color='green', label=r'$10x_c$')
plt.legend(title='Initial Position [m]')
plt.tight_layout()
plt.savefig('springftv.pdf')
plt.show()
