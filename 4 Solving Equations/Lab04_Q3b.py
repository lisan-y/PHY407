"""
PHY407 Lab 04
Question 3 b

Compares solutions and number of iterations to solve the nonlinear equation
x = 1 - np.exp(- c * x) using the methods of relaxation and overrelaxation

Author: Lisa Nasu-Yu, Sept 2021
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

# changes font of plots, but runs slower.
# mpl.rcParams['legend.frameon'] = False
# mpl.rcParams['figure.autolayout'] = True
#
# plt.rcParams.update({"text.usetex":True, "font.family": 'sans-serif', 'font.sans-serif': ['Helvetica']})
# plt.rcParams.update({"text.usetex":True, "font.family": 'serif', 'font.serif': ['Palatino']})


def f(c, x):
    return 1 - np.exp(- c * x)


def relaxation(f, c, accuracy):
    """
    Calculate solution to nonlinear equation using method of relaxation.
    :param f: [numpy function] equation to solve
    :param c: [float] constant
    :param accuracy: [float] minimum desired accuracy of solution
    :return: solution and number of iterations it took
    """
    error = 1.
    # choose initial x
    x1 = 1.
    iterations = 0
    while error > accuracy:
        x1, x2 = f(c, x1), x1
        error = np.abs(x2 - x1)
        iterations += 1
    return x1, iterations


def overrelaxation(f, c, accuracy, w):
    """
    Calculate solution to nonlinear equation using method of overelaxation.
    :param f: [numpy function] equation to solve
    :param c: [float] constant
    :param accuracy: [float] minimum desired accuracy of solution
    :param w: size of small increase in x
    :return: solution and number of iterations it took
    """
    error = 1.
    # choose initial x
    x1 = 1.
    iterations = 0
    while error > accuracy:
        x1, x2 = (1 + w) * f(c, x1) - w * x1, x1
        error = np.abs(x2 - x1)
        iterations += 1
    return x1, iterations


# compute solution and number of iterations to reach accuracy of 1e-6 with relaxation
x, it = relaxation(f, 2, 1e-6)
print('The solution to an accuracy of 1e-6 is {}. '
      '\nIt took {} iterations to reach this accuracy.'.format(x, it))

# compute solution and number of iterations to reach accuracy of 1e-6 with relaxation
w_array = np.arange(0, 1.1, 0.1)
x_array = []
it_array = []
for w in w_array:
    x, it = overrelaxation(f, 2, 1e-6, w)
    x_array.append(x)
    it_array.append(it)


np.savetxt("overrelaxaton.csv", np.array([w_array, x_array, it_array]).T, delimiter=' & ', newline=' \\\\\n')

