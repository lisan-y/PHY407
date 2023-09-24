"""
PHY407 Lab 03
Question 2

Numerically calculates the period of a particle on a spring, initially at rest.
The integral is computed with Gaussian quadrature.
The periods and errors are compared for various initial positions and order N.

Author: Lisa Nasu-Yu, Sept 2021
"""
import numpy as np
import matplotlib.pyplot as plt

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


# Define constants
c = 3e8  # speed of light [m/s]
k = 12  # spring constant [N/m]
m = 1  # mass [kg]


# Q2a

# Check classical limit of period for x0=1cm by comparing potential and kinetic energies
U = (1/2) * k * 0.01 ** 2
K = m * c ** 2

print("The ratio of potential energy and kinetic energy for a particle with x0=1 cm is U/K=", U/K)

# Calculate sampling points and period for N=8 and N=16 for x0=1cm
x, w = gaussxw(8)
T8, x8, w8 = period(0.01, x, w)

x, w = gaussxw(16)
T16, x16, w16 = period(0.01, x, w)

# Estimate fractional errors with e_N = I_2N - I_N
frac_error_8 = (T16 - T8) / T8

x, w = gaussxw(32)
T32, _, _ = period(0.01, x, w)
frac_error_16 = (T32 - T16) / T16

print("The period computed with Gaussian quadrature for N=8 is {} s with a fractional error of {}."
      .format(T8, frac_error_8))
print("The period computed with Gaussian quadrature for N=16 is {} s with a fractional error of {}."
      .format(T16, frac_error_16))

# Plot integrands at sampling points
plt.figure()
plt.title('Integrand Values')
plt.xlabel('Sampling points [cm]')
plt.ylabel(r'$4/g_k$ [s/m]')
plt.bar(x8 * 100, 4 / g(x8, 0.01), width=0.01, label='N=8')
plt.bar(x16 * 100, 4 / g(x16, 0.01), width=0.01, label='N=16')
plt.legend()
plt.tight_layout()
plt.savefig('2a_integrands.pdf')
plt.show()

# Plot weighted integrands at sampling points
plt.figure()
plt.title('Weighted Integrand Values')
plt.xlabel('Sampling points [cm]')
plt.ylabel(r'$4w_k/g_k$ [s/m]')
plt.bar(x8 * 100, 4 * w8 / g(x8, 0.01), width=0.01, label='N=8')
plt.bar(x16 * 100, 4 * w16 / g(x16, 0.01), width=0.01, label='N=16')
plt.legend()
plt.tight_layout()
plt.savefig('2a_weighted.pdf')
plt.show()


# Q2 b
# No code

# Q2 c
# Calculate sampling points and period for N=200 and N=400 for x0=1cm
x, w = gaussxw(200)
T200, _, _ = period(0.01, x, w)

x, w = gaussxw(400)
T400, _, _ = period(0.01, x, w)

# Estimate fractional errors for N=200
frac_error_200 = (T400 - T200) / T200 * 100

print("The period computed with Gaussian quadrature for N=200 is {} s with a percentage error of {}."
      .format(T200, frac_error_200))

# Plot T for 1 m < x0 < 10 * x_c, N=200
x0_range = np.logspace(np.log10(1), np.log10(10 * 3e8 / (2 * np.sqrt(3))))

plt.figure()
plt.title('Period for Particle on a Spring')
plt.xlabel(r'$x_0$ [m]')
plt.ylabel(r'T [s]')
plt.plot(x0_range, [period(i, x, w)[0] for i in x0_range])
plt.tight_layout()
plt.savefig('2c.pdf')
plt.show()
