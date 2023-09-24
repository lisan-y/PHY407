"""
PHY407 Lab 03
Question 3

Computes and plots the absolute errors of the derivative of the function
f(x)=exp(-x**2) for step size, h, ranging from 1e-16 to 1e0,
calculated by the forward and central difference schemes.

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


def f(x):
    """
    Computes the function exp(- x ** 2)
    :param x: [float] value at which the function is evaluated
    :return: [float] value of function at x
    """
    return np.exp(- x ** 2)


def forward_diff(x, h):
    """
    Computes the derivative of f(x)=exp(- x ** 2) using the forward difference scheme
    :param x: [float] value at which the function is evaluated
    :param h: [float] differentiation step size
    :return: [float] numerical value of derivative of f(x) at x
    """
    return (f(x + h) - f(x)) / h


def central_diff(x, h):
    """
    Computes the derivative of f(x)=exp(- x ** 2) using the central difference scheme
    :param x: [float] value at which the function is evaluated
    :param h: [float] differentiation step size
    :return: [float] numerical value of derivative of f(x) at x
    """
    return (f(x + h / 2) - f(x - h / 2)) / h


# Q 3a

# Compute numerical values for forward difference scheme
f_values = [forward_diff(0.5, h) for h in np.logspace(-16, 0, 17)]


# Q 3b

# Compute errors for forward difference scheme
f_errors = np.abs(-np.exp(-1 / 4) - f_values)
print('[\t\th\t\t\tnumerical value\t\terror\t\t]')
print(np.array([(np.logspace(-16, 0, 17)[i], f_values[i], f_errors[i]) for i in range(len(f_values))]).reshape((-1, 3)))

# np.savetxt("Q3ab.csv", np.array([(-16 + i, derivs[i], errors[i]) for i in range(len(derivs))]).reshape((-1, 3)), \
#            delimiter=' & ', fmt='%2.2e', newline=' \\\\\n')


# Q3 c, d

# Compute numerical values for central difference scheme
c_values = [central_diff(0.5, h) for h in np.logspace(-16, 0, 17)]
# Compute errors for central difference scheme
c_errors = np.abs(-np.exp(-1 / 4) - c_values)

# Plot errors
plt.figure()
plt.grid()
plt.title('Errors for Differentiation Schemes')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$h$ [log]')
plt.ylabel('Absolute Error [log]')
plt.plot(np.logspace(-16, 0, 17), f_errors, '^', label='Forward Difference')
plt.plot(np.logspace(-16, 0, 17), c_errors, 'v', label='Central Difference')
plt.legend()
plt.tight_layout()
# plt.savefig('3d.pdf')
plt.show()
