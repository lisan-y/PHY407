"""
PHY407 Lab 0
Question 2

Simulates time dependent Schrodinger Equation.

Author: Lisa Nasu-Yu, Nov 2021
"""
import numpy as np
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

# Define constants
h = 6.62607004e-34 / (2 * np.pi) # [m^2kg/s]
L = 1e-8  # [m]
m = 9.109e-31  # [kg]
sigma = L / 25
kappa = 500 / L

P = 1024   # x interval segments
xrange = np.linspace(-L / 2, L / 2, P)
dx = L / P
tau = 1e-18   # [s]
prange = np.arange(1, P, 1)

# Define initial condition of wavefunction (lab Eq. 21)
def psi_ic(x, x_0):
    # exponent part of psi(t=0)
    psi_exp = np.exp(-1 * (x - x_0)**2 / (4*sigma**2) + 1.j * kappa * x)

    # psi(t=0) is normalized so integral(conj(psi)*psi) = 1.
    # compute psi0 from this
    psi0 = np.sqrt(1 / (np.sum(psi_exp * np.conj(psi_exp)) * dx))
    return psi0 * psi_exp


# Define functions for potentials
def square_well(x):
    return 0


def harmonic_osc(x):
    w = 3e15
    return m * w**2 * x**2 / 2


def double_well(x):
    v0 = 6e-17  # [J]
    x1 = L / 4
    return v0 * ((x/x1)**2 - 1)**2


# Define discretized Hamiltonian matrix for a given potential
def H_D(V):
    A = -(h ** 2) / (2 * m * (dx ** 2))
    B = V(prange * dx - L / 2) - 2 * A
    matrix = np.eye(P - 1, k=0) * B
    matrix += np.eye(P - 1, k=1) * A
    matrix += np.eye(P - 1, k=-1) * A
    return matrix


# Define diagnostics functions as lab Eq. 16
def energy(psi, V):
    integrand = np.matmul(np.conjugate(psi), np.matmul(H_D(V), psi))
    # integrand = np.dot(np.conjugate(psi), np.dot(H_D(V), psi))
    return np.sum(integrand) * dx


def position(psi):
    integrand = np.sum(np.conjugate(psi) * xrange[1:] * psi)
    return np.sum(integrand) * dx


def normalization(psi):
    return np.sum(np.conjugate(psi) * psi) * dx


# Define L and R for Crank-Nicolson
def crank_nicolson(V):
    L = np.eye(P - 1, k=0) + (0.5 * 1.j * tau / h) * H_D(V)
    R = np.eye(P - 1, k=0) - (0.5 * 1.j * tau / h) * H_D(V)
    return L, R


# Define constants
x0 = L / 5
N = 3000  # time steps
potential = square_well

psi = psi_ic(xrange, x0)[1:]

L, R = crank_nicolson(potential)

# initial data arrays
energy_dat = [energy(psi, potential)]
x_dat = [position(psi)]
norm_dat = [normalization(psi)]
prob_density_dat = [np.conjugate(psi) * psi]

# time evolution
for i in range(N):
    # print progress
    if (i % (N//10)) == 0:
        print('Progress: ', i/N)
    # compute new psi
    psi = np.linalg.solve(L, np.matmul(R, psi))

    energy_dat.append(energy(psi, potential))
    x_dat.append(position(psi))
    norm_dat.append(normalization(psi))

    # store probability density if t=T/4 or T/2
    if i in [N//4, N//2]:
        prob_density_dat.append(np.conjugate(psi) * psi)
prob_density_dat.append(np.conjugate(psi) * psi)

# part a
plt.figure()
plt.title('Energy')
plt.xlabel('Time [s]')
plt.ylabel('Energy [J]')
plt.plot(np.arange(N+1) * tau, energy_dat)
plt.ylim(-2*energy_dat[1], 2*energy_dat[1])
plt.tight_layout()
# plt.savefig('energy_sw.pdf')
plt.show()

plt.figure()
plt.title(r'Normalization of $\psi$')
plt.xlabel('Time [s]')
plt.ylabel('Normalization')
plt.ylim(-2*norm_dat[1], 2*norm_dat[1])
plt.plot(np.arange(N+1) * tau, norm_dat)
plt.tight_layout()
# plt.savefig('norm_sw.pdf')
plt.show()

plt.figure()
plt.title('Expectation Position')
plt.xlabel('Time [s]')
plt.ylabel('Expected Position [m]')
plt.plot(np.arange(N+1) * tau, x_dat)
plt.tight_layout()
# plt.savefig('x_sw.pdf')
plt.show()

plt.figure()
plt.title('Probability Density')
plt.xlabel('Position [m]')
plt.ylabel('Probability Density [1/m]')
labels = [0, 'T/4', 'T/2', 'T']
print(np.shape(prob_density_dat))
for n in range(len(prob_density_dat)):
    plt.plot(xrange[1:], prob_density_dat[n], label=labels[n])
plt.legend(title='time')
plt.tight_layout()
# plt.savefig('prob_sw.pdf')
plt.show()
