import math as m
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

# mpl.rcParams['legend.frameon'] = False
# mpl.rcParams['figure.autolayout'] = True
#
# plt.rcParams.update({"text.usetex":True, "font.family": 'sans-serif', 'font.sans-serif': ['Helvetica']})
# plt.rcParams.update({"text.usetex":True, "font.family": 'serif', 'font.serif': ['Palatino']})

# Q3b
# 6.10 a


def f(c, x):
    return 1 - m.exp(- c * x)


def relaxation(f, c, accuracy):
    error = 1.
    # choose initial x
    x1 = 1.
    iterations = 0
    while error > accuracy:
        x1, x2 = f(c, x1), x1
        error = np.abs(x2 - x1)
        iterations += 1
    return x1, iterations


# 6.10 b
plt.figure()
plt.title(r'$x = 1 - e^{-cx}$')
plt.xlabel('c')
plt.ylabel('x')
for i in np.arange(0, 3, 0.01):
    plt.plot(i, relaxation(f, i, 1e-6)[0], 'b.')
plt.tight_layout()
# plt.savefig('relaxation.pdf')
plt.show()


x, it = relaxation(f, 2, 1e-6)
print('The solution to an accuracy of 1e-6 is {}. '
      '\nIt took {} iterations to reach this accuracy.'.format(x, it))


