"""
PHY407 Lab 07
Question 2

Compute energies of first 3 states of harmonic quantum oscillator

Author: Lisa Nasu-Yu, Oct 2021
"""
from numpy import array, arange, abs, exp, pi
import matplotlib.pyplot as plt

# Constants
m = 9.1094e-31     # Mass of electron
hbar = 1.0546e-34  # Planck's constant over 2*pi
e = 1.6022e-19     # Electron charge
L = 5.2918e-11     # Bohr radius
V0 = 50 * e  # [J]
a = 1e-11  # [m]
x0 = -10 * a
xf = 10 * a
N = 1000
h = (xf - x0) / N

# Potential functions
def V2(x):
    return V0 * x**2 / a**2


def V4(x):
    return V0 * x ** 4 / a ** 4


def f(r,x,E, V):
    # adapted from Newman
    psi = r[0]
    phi = r[1]
    fpsi = phi
    fphi = (2*m/hbar**2)*(V(x)-E)*psi
    return array([fpsi,fphi],float)


# Calculate the wavefunction for a particular energy
def solve(E, V):
    # adapted from Newman
    psi = 0.0
    phi = 1.0
    r = array([psi, phi], float)
    psi_arr = []

    # RK4
    for x in arange(x0, xf, h):
        psi_arr.append(r[0])
        k1 = h*f(r, x, E, V)
        k2 = h*f(r+0.5*k1, x+0.5*h, E, V)
        k3 = h*f(r+0.5*k2, x+0.5*h, E, V)
        k4 = h*f(r+k3, x+h, E, V)
        r += (k1+2*k2+2*k3+k4)/6

    return psi_arr


# Main program to find the energy using the secant method
def secant(E1, E2, V):
    # adapted from Newman
    psi_arr2 = solve(E1, V)
    psi2 = psi_arr2[-1]
    target = e / 1000
    while abs(E1 - E2) > target:
        psi1 = psi2
        psi_arr2 = solve(E2, V)
        psi2 = psi_arr2[-1]
        E1, E2 = E2, E2 - psi2 * (E2 - E1) / (psi2 - psi1)

    print("E =", E2 / e, "eV")
    return E2/e, psi_arr2


def trapezoid(psi):
    """
    Evalutes probability density using trapezoid rule with N slices for an array of points of psi
    :param f: [array] integrand
    :param a: [float] lower bound of integration
    :param b: [float] upper bound of integration
    :param N: [int] number of slices
    :return: [float] solution to integral with trapezoid rule for N slices
    """
    # integrate left half
    s = 0.5 * abs(psi[0])**2 + 0.5 * abs(psi[N // 2])**2  # end slices
    for k in range(1, N // 2):  # add each interior slice
        s += abs(psi[k])**2
    # multiply by 2 for full unnormalized probability density
    prob_density = 2 * h * s
    return prob_density**0.5


# 3 energy states for V(x) as given in a)
E0a, psi0a = secant(0, e, V2)
E1a, psi1a = secant(200*e, 400*e, V2)
E2a, psi2a = secant(500*e, 700*e, V2)

# check that energies are evenly separated
print('Change in spacing of energy levels: {}eV'.format((E0a - E1a) - (E1a - E2a)))

# b)
# 3 energy states for V(x) as given in b)
E0b, psi0b = secant(0, e, V4)
E1b, psi1b = secant(400*e, 600*e, V4)
E2b, psi2b = secant(900*e, 1100*e, V4)

# c)
# find normalization factors
a0 = trapezoid(psi0a)
a1 = trapezoid(psi1a)
a2 = trapezoid(psi2a)

b0 = trapezoid(psi0b)
b1 = trapezoid(psi1b)
b2 = trapezoid(psi2b)

# analytic wavefunctions
xpoints = arange(x0, xf, h)


# functions to compute probability densities for 1st 3 states
def psi0_analytic(x, alpha):
    y = x * alpha**0.5
    return (alpha/pi)**0.25 * exp(-y**2 / 2)


def psi1_analytic(x, alpha):
    y = x * alpha**0.5
    return (alpha/pi)**0.25 * 2**0.5 * y * exp(-y**2 / 2)


def psi2_analytic(x, alpha):
    y = x * alpha**0.5
    return (alpha/pi)**0.25 * 2**(-0.5) * (2*y**2 - 1) * exp(-y**2 / 2)


# harmonic
# 1st level
w = 2 * E0a*e / hbar
alpha = m * w / hbar
analytic0a = psi0_analytic(xpoints, alpha)

# 2nd level
w = (2/3) * E1a*e / hbar
alpha = m * w / hbar
analytic1a = psi1_analytic(xpoints, alpha)

# 3rd level
w = (2/5) * E2a*e / hbar
alpha = m * w / hbar
analytic2a = psi2_analytic(xpoints, alpha)

# plot for x in range -5a to 5a
startx = N // 4
endx = 3 * N // 4
plt.figure()
plt.title('Harmonic Oscillator', fontsize=16)
plt.plot(xpoints[startx:endx], psi0a[startx:endx] / a0, color='blue', label='0')
plt.plot(xpoints[startx:endx], analytic0a[startx:endx], '--', color='cyan', label='0 (analytic)')
plt.plot(xpoints[startx:endx], psi1a[startx:endx] / -a1, color='orange', label='1')
plt.plot(xpoints[startx:endx], analytic1a[startx:endx], '--', color='yellow', label='1 (analytic)')
plt.plot(xpoints[startx:endx], psi2a[startx:endx] / a2, color='green', label='2')
plt.plot(xpoints[startx:endx], analytic2a[startx:endx], '--', color='lime', label='2 (analytic)')
plt.xlabel('x (m)', fontsize=14)
plt.ylabel(r'$\psi$', fontsize=14)
plt.legend(title='Energy Level', fontsize=14)
plt.tight_layout()
plt.savefig('harmonic.pdf')
plt.show()


# anharmonic
plt.figure()
plt.title('Anharmonic Oscillator', fontsize=16)
plt.plot(xpoints[startx:endx], psi0b[startx:endx] / b0, color='blue', label='0')
plt.plot(xpoints[startx:endx], psi1b[startx:endx] / -b1, color='orange', label='1')
plt.plot(xpoints[startx:endx], psi2b[startx:endx] / b2, color='green', label='2')
plt.xlabel('x (m)', fontsize=14)
plt.ylabel(r'$\psi$', fontsize=14)
plt.legend(title='Energy Level', fontsize=14)
plt.tight_layout()
plt.savefig('anharmonic.pdf')
plt.show()
