"""
PHY407 Lab 1 Q1 d

Computes and plots the position (x and y) and x and y components of the
velocity (vx and vy) of the planet Mercury as a function of time under
relativistic gravitational force given the initial conditions for x, y, vx and vy and
using the Euler-Cromer method.
Author: Lisa Nasu-Yu, Sept. 2021
"""
import numpy as np
import matplotlib.pyplot as plt


# Define constants
M_s = 1  # mass of Sun [solar mass]
G = 39.5  # gravitational constant [AU^3 * M_s^-1 * yr^-2]
dt = 0.0001  # time step [yr]
t_f = 3  # duration of integration [yr]

# Initiate time, position, and velocity arrays
t = np.arange(0, t_f, dt)  # time [yr]
x = np.zeros(len(t))
vx = np.zeros(len(t))
y = np.zeros(len(t))
vy = np.zeros(len(t))

# Initial conditions:
x[0] = 0.47  # [AU]
y[0] = 0  # [AU]
vx[0] = 0  # [AU/yr]
vy[0] = 8.17  # [AU/yr]

# Calculate position & velocity by integrating
for i in range(len(t) - 1):
    a = 0.01  # constant [AU**2]
    r = np.sqrt(x[i] ** 2 + y[i] ** 2)  # distance between Mercury and Sun [AU]
    grav_eqn = (G * M_s / r ** 3) * (1 + a / r ** 2)
    vx[i + 1] = vx[i] - grav_eqn * x[i] * dt
    vy[i + 1] = vy[i] - grav_eqn * y[i] * dt
    x[i + 1] = x[i] + vx[i + 1] * dt
    y[i + 1] = y[i] + vy[i + 1] * dt

# plot x vs y position
plt.figure()
plt.rc('font', size=14)
plt.axis('equal')
plt.ylabel('Y Position [AU]')
plt.xlabel('X Position [AU]')
plt.title('X vs Y Positions of Mercury \nwith General Relativity Orbits')
plt.plot(x, y)
plt.savefig('Q1d.pdf')
plt.show()

plt.figure()
plt.rc('font', size=14)
plt.title('Velocity Components of Mercury \nwith General Relativity')
plt.ylabel('Velocity [AU/yr]')
plt.xlabel('Time [yr]')
plt.plot(t, vx, label=r'$v_x$')
plt.plot(t, vy, label=r'$v_y$')
plt.legend()
plt.savefig('Q1d_v.pdf')
plt.show()
