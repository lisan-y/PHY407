"""
PHY407 Lab 1 Q1 b & c

Q1 b
    Pseudocode for computing and plotting the position (x and y) and x and y components of the
    velocity (vx and vy) of a planet as a function of time under
    Newtonian gravity force given the initial conditions for x, y, vx and vy and
    using the Euler-Cromer method.
    Author: Emma Jarvis, Sept. 2021

Q1 c
    Codes the 1b pseudocode for the planet Mercury, a duration of 1 year, and a time step of
    0.0001 year.
    Computes and plots the angular momentum of Mercury using the calculated velocities and positions.
    Author: Lisa Nasu-Yu, Sept. 2021
"""
import numpy as np
import matplotlib.pyplot as plt

# Q1 b

# Pseudocode:
# Import necessary packages (numpy and matplotlib.pyplot)
# 1. Define the gravitational constant, G and the mass of the sun, M_s.
# 2. Define initial x0, y0, vx0, vy0.
# 3. Define time step dt, and final time t_f.
# 4. Initialize the time array from initial time to final time with a time step of dt.
# 5. Initialize x, y, vx, and vy arrays with the same number of elements as the time array.
# 6. Set first element of x, y, vx, and vy with their respective initial values.
# 7. For number of iterations corresponding to 1 less than the length of the time array:
# Increment the vx array with vx[i+1] = vx[i] - (G*M_s/sqrt(x[i]**2 + y[i]**2)**3) * x[i] * dt
# Increment the vy array with vy[i+1] = vy[i] - (G*M_s/sqrt(x[i]**2 + y[i]**2)**3) * y[i] * dt
# Increment the x array with Euler-Cromer: x[i+1] = x[i] + vx[i+1] * dt
# Increment the y array with Euler-Cromer: y[i+1] = y[i] + vy[y+1] * dt
# 7. Plot vx vs. t
# 8. Plot vy vs. t
# 9. Plot x vs. y

##########################################################

# Q1 c

# Define constants
M_s = 1  # mass of Sun [solar mass]
G = 39.5  # gravitational constant [AU^3 * M_s^-1 * yr^-2]
dt = 0.0001  # time step [yr]
t_f = 1  # duration of integration [yr]

# Initiate time, position, and velocity arrays
t = np.arange(0, t_f, dt)  # time [yr]
x = np.zeros(len(t))
vx = np.zeros(len(t))
y = np.zeros(len(t))
vy = np.zeros(len(t))

# Set initial conditions in arrays
x[0] = 0.47  # [AU]
y[0] = 0  # [AU]
vx[0] = 0  # [AU/yr]
vy[0] = 8.17  # [AU/yr]

# Calculate position & velocity
for i in range(len(t) - 1):
    grav_eqn = (G * M_s / np.sqrt(x[i] ** 2 + y[i] ** 2) ** 3)
    vx[i + 1] = vx[i] - grav_eqn * x[i] * dt
    vy[i + 1] = vy[i] - grav_eqn * y[i] * dt
    x[i + 1] = x[i] + vx[i + 1] * dt
    y[i + 1] = y[i] + vy[i + 1] * dt

# plot time vs x and y velocities
plt.figure()
plt.rc('font', size=14)
plt.title('Velocity Components of Mercury \nunder Newtonian Orbits')
plt.ylabel('Velocity [AU/yr]')
plt.xlabel('Time [yr]')
plt.plot(t, vx, label=r'$v_x$')
plt.plot(t, vy, label=r'$v_y$')
plt.legend()
plt.savefig('Q1c_v.pdf')
plt.show()

# plot x vs y position
plt.figure()
plt.rc('font', size=14)
plt.axis('equal')
plt.ylabel('Y Position [AU]')
plt.xlabel('X Position [AU]')
plt.title('X vs Y Positions of Mercury \nunder Newtonian Orbits')
plt.plot(x, y)
plt.savefig('Q1c_xy.pdf')
plt.show()

# Calculate angular momentum
M_m = 1.651  # mass of Mercury [solar mass]
# Angular momentum, L = m*v*r [AU * M_S ** 2 / yr]
L = M_m * np.sqrt(vx ** 2 + vy ** 2) * np.sqrt(x ** 2 + y ** 2)

# plot t vs angular momentum
plt.figure()
plt.rc('font', size=11)
plt.title('Angular Momentum of Mercury under Newtonian Orbits')
plt.ylabel(r'Angular Momentum [$M_sAU^2yr^{-1}]$')
plt.xlabel('Time [yr]')
plt.plot(t, L)
plt.savefig('Q1c_L.pdf')
plt.show()

print('Q1c\nMax x velocity: {:.03g} AU/yr'
      '\nMin x velocity of Mercury: {:.03g} AU/yr'
      '\nx velocity Amplitude {:.03g} AU/yr'.format(max(vx), min(vx), (max(vx) - min(vx))/2))

print('Q1c\nMax y velocity: {:.03g} AU/yr'
      '\nMin y velocity of Mercury: {:.03g} AU/yr'
      '\ny velocity Amplitude {:.03g} AU/yr'.format(max(vy), min(vy), (max(vy) - min(vy))/2))

print('Q1c\nMax Angular Momentum of Mercury: {:.03g} AU*M_s/yr'
      '\nMin Angular Momentum of Mercury: {:.03g} AU*M_s/yr'.format(max(L), min(L)))
