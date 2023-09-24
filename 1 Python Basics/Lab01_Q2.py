"""
PHY407 Lab 1 Q2

Given the initial conditions for x, y, vx, and vy, of objects in a 3 body system centred
at the Sun, compute the position (x and y) and x and y components of the velocity (vx and vy)
of the remaining 2 objects using the Euler-Cromer method.
Plot the spatial coordinates in the x-y planes of the 2 objects over some period of time.

Q2 a
    Simulate orbits of Earth and Jupiter over 10 years
Q2 b
    Simulate orbits of Earth and a Jupiter-like object with a mass of 1 solar mass,  over 3 years
Q2 c
    Simulate orbits of asteroid and Jupiter over 20 years
Author: Lisa Nasu-Yu, Sept. 2021
"""

import numpy as np
import matplotlib.pyplot as plt

# Q2 a

# Pseudocode:
# Import necessary packages (numpy and matplotlib.pyplot)
# 1. Define the gravitational constant, G and the mass of the sun, M_s

# 2. Apply Q1 b pseudocode to calculate x & y position and velocity for Jupiter,
#   and corresponding initial conditions xj0, yj0, vxj0, vyj0

# 3. Define time step, dt=0.0001 year, and duration of integration, t_f=10 years
# 4. Define Earth's initial x & y position and velocity: xe0, ye0, vxe0, vye0
# 5. Initialize time array from initial time to final time with a time step of dt.
# 6. Initialize x & y position and velocity (x, y, vx, and vy) arrays with the same
#   number of elements as the time array.
# 7. Set first element of x, y, vx, and vy with their respective initial values.
# 8. For number of iterations corresponding to 1 less than the length of the time array:
#         Calculate distance between Earth and Jupiter, rj = sqrt((xj[i] - x[i])**2 + (yj[i] - y[i])**2)
#         Define gravitational force fj on Earth due to Jupiter (G*M_s/rj**3)
#         Calculate distance between Earth and Sun, rs = sqrt(x[i]**2 + y[i]**2)
#         Define gravitational force fj on Earth due to Jupiter (G*M_s/rs**3)
#         Total gravitational force felt by Earth is the sum of that from Jupiter and the Sun:
#           Increment the vx array with vx[i+1] = vx[i] - (fs * x[i] - fj * (xj[i] - x[i])) * dt
#           Increment the vy array with vy[i+1] = vy[i] - (fs * y[i] - fj * (yj[i] - y[i])) * dt
#           Increment the x array with Euler-Cromer: x[i+1] = x[i] + vx[i+1] * dt
#           Increment the y array with Euler-Cromer: y[i+1] = y[i] + vy[y+1] * dt
# 9. Plot x vs. y

########################################################################################

# Functions


def initiate_arrays(x0, y0, vx0, vy0, dt, t_f):
    """
    Creates an array for time as numpy array from 0 to t, with step of size dt.
    Initiates an array of zeroes of length t for x, y, vx, and vy.

    :param x0: float, initial x position [AU]
    :param y0: float, initial y position [AU]
    :param vx0: float, initial x velocity [AU/yr]
    :param vy0: float, initial y velocity [AU/yr]
    :param dt: float, timestep [yr]
    :param t_f: float or int, final time [yr]
    :return: arrays for x position, y position, x velocity, y velocity
    """

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
    """
    Given the the ith x and y positions and velocities,
    compute the the (i + 1)th x and y positions and velocities
    using the the Euler-Cromer method.

    :param grav_eqn: numpy equation, gravitational force equation * r
    :param xi: float, ith x position
    :param yi: float, ith y position
    :param vxi: float, ith x velocity
    :param vyi: float, ith y velocity
    :param dt: float, time step [yr]
    :return: (i + 1)th x and y positions and velocities
    """
    vx_iplus1 = vxi - grav_eqn * xi * dt
    vy_iplus1 = vyi - grav_eqn * yi * dt
    x_iplus1 = xi + vx_iplus1 * dt
    y_iplus1 = yi + vy_iplus1 * dt
    return x_iplus1, y_iplus1, vx_iplus1, vy_iplus1


def three_body(dt, t_f, x0, y0, vx0, vy0, xj0, yj0, vxj0, vyj0, M_j):
    """
    Simulate orbits of 2 objects in a 3 body system
    centred at the sun using the Euler-Cromer method.

    :param dt: float, time step [yr]
    :param t_f: float, final time [yr]
    :param x0: float, initial x position of object 1 [AU]
    :param y0: float, initial y position of object 1 [AU]
    :param vx0: float, initial x velocity of object 1 [AU/yr]
    :param vy0: float, initial y velocity of object 1 [AU/yr]
    :param xj0: float, initial x position of object 2 [AU]
    :param yj0: float, initial y position of object 2 [AU]
    :param vxj0: float, initial x velocity of object 2 [AU/yr]
    :param vyj0: float, initial y velocity of object 2 [AU/yr]
    :param M_j: Mass [solar mass]
    :return: arrays for x and y position object 1 and 2
    """
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

    return x, y, xj, yj


def plot_xy(x1, x2, y1, y2, obj1, obj2, filename, time=''):
    plt.figure()
    plt.rc('font', size=14)
    plt.axis('equal')
    plt.ylabel('Y Position [AU]')
    plt.xlabel('X Position [AU]')
    plt.title(obj1 + ' and ' + obj2 + ' Orbits ' + time)
    plt.plot(x1, y1, label=obj1)
    plt.plot(x2, y2, label=obj2)
    plt.legend()
    plt.savefig(filename)
    plt.show()

###############################################################################################

# Q2a code (pseudocode above)

# Define constants
M_s = 1  # mass of Sum [solar mass]
M_j = 1e-3 * M_s  # mass of Jupiter [solar mass]
G = 39.5  # gravitational constant [AU^3 * M_s^-1 * yr^-2]
dt = 0.0001  # time step [yr]
t_f = 10  # duration of integration [yr]

# Set Jupiter's initial conditions:
xj0 = 5.2  # [AU]
yj0 = 0  # [AU]
vxj0 = 0  # [AU/yr]
vyj0 = 2.63  # [AU/yr]

# Set Earth's initial conditions:
x0 = 1.0  # [AU]
y0 = 0  # [AU]
vx0 = 0  # [AU/yr]
vy0 = 6.18  # [AU/yr]


# compute spatial coordinates in x-y plane for Jupiter and Earth
x_a, y_a, xj_a, yj_a = three_body(dt, t_f, x0, y0, vx0, vy0, xj0, yj0, vxj0, vyj0, M_j)
# plot
plot_xy(x_a, xj_a, y_a, yj_a, 'Earth', 'Jupiter', 'Q2_a.pdf')


# Q2 b

# Change Jupiter mass to 1 solar mass
M_j = M_s
# Change duration to 3 years
t_f = 3  # [yr]

# compute spatial coordinates in x-y plane for high mass Jupiter and Earth
x_b1, y_b1, xj_b1, yj_b1 = three_body(dt, t_f, x0, y0, vx0, vy0, xj0, yj0, vxj0, vyj0, M_j)
# plot
plot_xy(x_b1, xj_b1, y_b1, yj_b1, 'Earth', 'High Mass Jupiter', 'Q2_b_3y.pdf', time='3 Years')


# Change duration to 6 years
t_f = 6  # [yr]

# compute spatial coordinates in x-y plane for high mass Jupiter and Earth
x_b2, y_b2, xj_b2, yj_b2 = three_body(dt, t_f, x0, y0, vx0, vy0, xj0, yj0, vxj0, vyj0, M_j)
# plot
plot_xy(x_b2, xj_b2, y_b2, yj_b2, 'Earth', 'High Mass Jupiter', 'Q2_b_6y.pdf', time='6 Years')


# Q2 c

# Change mass of Jupiter to 1e-3 solar mass
M_j = 1e-3 * M_s

# Set asteroid's initial conditions:
xa0 = 3.3  # [AU]
ya0 = 0  # [AU]
vax0 = 0  # [AU/yr]
vay0 = 3.46  # [AU/yr]

# set duration to 20 years
t_f = 20  # [yrs]


# compute spatial coordinates in x-y plane for asteroid and Jupiter
x_c, y_c, xj_c, yj_c = three_body(dt, t_f, xa0, ya0, vax0, vay0, xj0, yj0, vxj0, vyj0, M_j)
# plot
plot_xy(x_c, xj_c, y_c, yj_c, 'Asteroid', 'Jupiter', 'Q2_c.pdf')
