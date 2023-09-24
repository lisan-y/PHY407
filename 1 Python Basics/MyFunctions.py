'''
PHY407
Module containing all functions used in Lab 1
'''
import numpy as np


M_s = 1  # mass of Sun [solar mass]
G = 39.5  # gravitational constant [AU^3 * M_s^-1 * yr^-2]


def initiate_arrays(x0, y0, vx0, vy0, dt, t_f):

    t = np.arange(0, t_f, dt)  # time array with time step dt [yr]

    # Initialize x & y position and velocity arrays of length t
    x = np.zeros(len(t))
    vx = np.zeros(len(t))
    y = np.zeros(len(t))
    vy = np.zeros(len(t))

    # Set initial conditions in arrays
    vx[0] = vx0
    x[0] = x0
    vy[0] = vy0
    y[0] = y0

    return t, x, y, vx, vy


def increment(grav_eqn, xi, yi, vxi, vyi, dt):
    vx_iplus1 = vxi - grav_eqn * xi * dt
    vy_iplus1 = vyi - grav_eqn * yi * dt
    x_iplus1 = xi + vx_iplus1 * dt
    y_iplus1 = yi + vy_iplus1 * dt
    return x_iplus1, y_iplus1, vx_iplus1, vy_iplus1

def newtonian_orbit(x0, y0, vx0, vy0, dt, t_f):

    # Initial conditions:
    x0 = x0  # [AU]
    y0 = y0  # [AU]
    vx0 = vx0  # [AU/yr]
    vy0 = vy0  # [AU/yr]

    # Set time step and time array
    dt = dt  # [yr]
    t = np.arange(0, t_f, dt)  # [yr]

    # Initialize x & y position and velocity arrays of length t
    x = np.zeros(len(t))
    vx = np.zeros(len(t))
    y = np.zeros(len(t))
    vy = np.zeros(len(t))

    # Set initial conditions in arrays
    vx[0] = vx0
    x[0] = x0
    vy[0] = vy0
    y[0] = y0

    # Calculate position & velocity by integrating
    for i in range(len(t) - 1):
        grav_eqn = (G * M_s / np.sqrt(x[i] ** 2 + y[i] ** 2) ** 3)
        vx[i + 1] = vx[i] - grav_eqn * x[i] * dt
        vy[i + 1] = vy[i] - grav_eqn * y[i] * dt
        x[i + 1] = x[i] + vx[i + 1] * dt
        y[i + 1] = y[i] + vy[i + 1] * dt

    return x, y, vx, vy, t


def three_body(dt, t_f, x0, y0, vx0, vy0, xj0, yj0, vxj0, vyj0, M_j):
    # Simulate Jupiter's orbit
    # initiate time, position, and velocity arrays for Jupiter
    t, xj, yj, vxj, vyj = initiate_arrays(xj0, yj0, vxj0, vyj0, dt, t_f)

    # Calculate Jupiter's position & velocity by integrating
    for i in range(len(t) - 1):
        grav_eqn = (G * M_s / np.sqrt(xj[i] ** 2 + yj[i] ** 2) ** 3)
        xj[i+1], yj[i+1], vxj[i+1], vyj[i+1] = increment(grav_eqn, xj[i], yj[i], vxj[i], vyj[i], dt)

    # Simulate Earth's orbit
    # initiate time, position, and velocity arrays for Earth
    t, x, y, vx, vy = initiate_arrays(x0, y0, vx0, vy0, dt, t_f)

    # Calculate Earth's position & velocity, summing force of gravity from Jupiter and Sun
    for i in range(len(t) - 1):
        rj = np.sqrt((xj[i] - x[i])**2 + (yj[i] - y[i])**2)  # distance between Jupiter and Earth [AU]
        rs = np.sqrt(x[i] ** 2 + y[i] ** 2)  # distance between Sun and Earth [AU]
        fj = G*M_j/rj**3  # component of gravitational force from Jupiter
        fs = G*M_s/rs**3  # component of gravitational force from Sun
        vx[i + 1] = vx[i] - (fs * x[i] + fj * (xj[i] - x[i])) * dt
        vy[i + 1] = vy[i] - (fs * y[i] + fj * (yj[i] - y[i])) * dt
        x[i + 1] = x[i] + vx[i + 1] * dt
        y[i + 1] = y[i] + vy[i + 1] * dt

    return x, y, vx, vy, xj, yj, vxj, vyj
