"""
PHY407 Lab 07
Question 1

Use adaptive step RK4 to simulate trajectory of a ball around a space rod.

Author: Lisa Nasu-Yu, Oct 2021
"""
# Adapted from solution to Newman 8.8, Space garbage, authored by Nico Grisouard, Univ. of Toronto

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import time


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

# plot settings
ftsz = 16
font = {'family': 'normal', 'size': ftsz}  # font size
rc('font', **font)

# %% This next part adapted from Newman's odesim.py --------------------------|
a = 0.0
b = 10.0

# Q1 a & b time adaptive step rk4
h = 0.01
delta = 1e-6  # target error per second [m/s]

adaptive_tpoints = []
adaptive_xpoints = []
adaptive_vxpoints = []  # the future dx/dt
adaptive_ypoints = []
adaptive_vypoints = []  # the future dy/dt

# time adaptive step rk4
# begin timer
timer_start = time.time()

# below: ordering is x, dx/dt, y, dy/dt
r = np.array([1., 0., 0., 1.], float)
t = a

while t < b:
    rho = 0
    adaptive_tpoints.append(t)
    adaptive_xpoints.append(r[0])
    adaptive_vxpoints.append(r[1])
    adaptive_ypoints.append(r[2])
    adaptive_vypoints.append(r[3])
    while rho < 1.:
        # step 1 with h
        k1 = h*rhs(r)  # all the k's are vectors
        k2 = h*rhs(r + 0.5*k1)  # note: no explicit dependence on time of the RHSs
        k3 = h*rhs(r + 0.5*k2)
        k4 = h*rhs(r + k3)
        rh1 = r + (k1 + 2*k2 + 2*k3 + k4)/6
        # step 2 with h
        k1 = h*rhs(rh1)  # all the k's are vectors
        k2 = h*rhs(rh1 + 0.5*k1)  # note: no explicit dependence on time of the RHSs
        k3 = h*rhs(rh1 + 0.5*k2)
        k4 = h*rhs(rh1 + k3)
        rh2 = rh1 + (k1 + 2*k2 + 2*k3 + k4)/6
        # repeat with step size 2*h
        k1 = 2*h*rhs(r)  # all the k's are vectors
        k2 = 2*h*rhs(r + 0.5*k1)  # note: no explicit dependence on time of the RHSs
        k3 = 2*h*rhs(r + 0.5*k2)
        k4 = 2*h*rhs(r + k3)
        r2h = r + (k1 + 2*k2 + 2*k3 + k4)/6
        # compute ratio
        ex = (rh2[0] - r2h[0]) / 30
        ey = (rh2[2] - r2h[2]) / 30
        rho = h * delta / np.sqrt(ex**2 + ey**2)
        if rho > 1.:
            t += 2 * h
        # change step size by rho**1/4 or limiting factor 2
        h *= min(1.5, rho**0.25)
    r = rh2

# end timer
adaptive_timer = time.time() - timer_start
print('Adaptive step RK4 time taken: {}s'.format(adaptive_timer))

# Q1 b time uniform step rk4
N = 10000
h = (b-a)/N

tpoints = np.arange(a, b, h)
xpoints = []
vxpoints = []  # the future dx/dt
ypoints = []
vypoints = []  # the future dy/dt

# begin timer
timer_start = time.time()
# below: ordering is x, dx/dt, y, dy/dt
r = np.array([1., 0., 0., 1.], float)
for t in tpoints:
    xpoints.append(r[0])
    vxpoints.append(r[1])
    ypoints.append(r[2])
    vypoints.append(r[3])
    k1 = h * rhs(r)  # all the k's are vectors
    k2 = h * rhs(r + 0.5 * k1)  # note: no explicit dependence on time of the RHSs
    k3 = h * rhs(r + 0.5 * k2)
    k4 = h * rhs(r + k3)
    r += (k1 + 2 * k2 + 2 * k3 + k4) / 6

# end timer
uniform_timer = time.time() - timer_start
print('Uniform step RK4 time taken: {}s'.format(uniform_timer))

# plot both methods
plt.figure()
plt.plot(xpoints, ypoints, ':', label='Uniform step')
plt.plot(adaptive_xpoints, adaptive_ypoints, 'k.', label='Adaptive step')
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title('Trajectory of Ball Around a Space Rod', fontsize=ftsz)
plt.legend()
plt.tight_layout()
plt.savefig('Garbage.pdf')
plt.show()

# Q1c plot timestep
dtpoints = np.array(adaptive_tpoints[1:]) - np.array(adaptive_tpoints[:-1])
plt.figure()
plt.title('Time Step Size')
plt.xlabel('Time [s]')
plt.ylabel('Step Size [s]')
plt.plot(adaptive_tpoints[3:-1], dtpoints[3:])
plt.tight_layout()
plt.savefig('timesteps.pdf')
plt.show()



