# Solution to Newman 8.8, Space garbage.
# Author: Nico Grisouard, Univ. of Toronto

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


def rhs(r):
    """ The right-hand-side of the equations
    INPUT:
    r = [x, vx, y, vy] are floats (not arrays)
    note: no explicit dependence on time
    OUTPUT:
    1x2 numpy array, rhs[0] is for x, rhs[1] is for vx, etc"""
    M = 10.
    L = 2.

    x = r[0]
    vx = r[1]
    y = r[2]
    vy = r[3]

    r2 = x**2 + y**2
    Fx, Fy = - M * np.array([x, y], float) / (r2 * np.sqrt(r2 + .25*L**2))
    return np.array([vx, Fx, vy, Fy], float)


ftsz = 16
font = {'family': 'normal', 'size': ftsz}  # font size
rc('font', **font)

# %% This next part adapted from Newman's odesim.py --------------------------|
a = 0.0
b = 10.0
N = 1000  # let's leave it at that for now
h = (b-a)/N
delta = 1e-6

tpoints = []
xpoints = []
vxpoints = []  # the future dx/dt
ypoints = []
vypoints = []  # the future dy/dt

# below: ordering is x, dx/dt, y, dy/dt
r = np.array([1., 0., 0., 1.], float)
t = 0.
rho = 0
while t < b:
    xpoints.append(r[0])
    vxpoints.append(r[1])
    ypoints.append(r[2])
    vypoints.append(r[3])
    while rho < 1.:
        k1 = h*rhs(r)  # all the k's are vectors
        k2 = h*rhs(r + 0.5*k1)  # note: no explicit dependence on time of the RHSs
        k3 = h*rhs(r + 0.5*k2)
        k4 = h*rhs(r + k3)
        rh1 = r + (k1 + 2*k2 + 2*k3 + k4)/6

        k1 = h*rhs(rh1)  # all the k's are vectors
        k2 = h*rhs(rh1 + 0.5*k1)  # note: no explicit dependence on time of the RHSs
        k3 = h*rhs(rh1 + 0.5*k2)
        k4 = h*rhs(rh1 + k3)
        rh2 = rh1 + (k1 + 2*k2 + 2*k3 + k4)/6

        k1 = 2*h*rhs(r)  # all the k's are vectors
        k2 = 2*h*rhs(r + 0.5*k1)  # note: no explicit dependence on time of the RHSs
        k3 = 2*h*rhs(r + 0.5*k2)
        k4 = 2*h*rhs(r + k3)
        r2h = r + (k1 + 2*k2 + 2*k3 + k4)/6

        ex = (rh2[0] - r2h[0]) / 30
        ey = (rh2[2] - r2h[2]) / 30
        print(ex, ey)
        rho = h * delta / np.sqrt(ex**2 + ey**2)
        h *= rho**0.25
    t += 2 * h
    r = rh2

plt.figure()
plt.plot(xpoints, ypoints, ':')
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title('Trajectory of a ball bearing around a space rod.', fontsize=ftsz)
plt.axis('equal')
plt.grid()
plt.tight_layout()
# plt.savefig('Garbage.png', dpi=150)
plt.show()
