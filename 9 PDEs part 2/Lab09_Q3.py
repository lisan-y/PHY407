""" Computing the vibrations of a string, hit by a hammer
Author: Nico Grisouard, University of Toronto
Inspired by Newman's heat.py """
##with section added to include the LW method


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

ftsz = 9
font = {'family': 'DejaVu Sans', 'size': ftsz}  # font size
rc('font', **font)
dpi = 100

# %% Constants and adjustable parameters -------------------------------------|
Lx = 2.*np.pi  # [m] length of domain
epsilon = 1.  # [] the non-linear parameter
Dt = 5e-3  # [s] time step
Dx = 2e-2  # [m] Grid step
Tf = 2.  # [s] duration of integration

t_record = [0.25, 0.5, 1., 1.5]  # the times of the snapshots
delta = Dt*0.5  # tolerance to detect when t_record is hit

# %% Dependent parameters ----------------------------------------------------|
beta = epsilon*Dt/Dx

# Independent variables
x = np.arange(0, Lx+Dx, Dx)
Nx = len(x)

t = np.arange(0, Tf+Dt, Dt)
Nt = len(t)

# Create arrays
um = np.sin(x)  # mm for minus, i.e. j-1
u = np.empty(Nx, float)  # time step j
up = np.empty(Nx, float)  # p for plus, time step j+1
um[0] = 0.  # boundary conditions
um[-1] = 0.  # boundary conditions
u[0] = 0.  # boundary conditions
u[-1] = 0.  # boundary conditions
up[0] = 0.  # boundary conditions
up[-1] = 0.  # boundary conditions

u_record = np.empty((len(t_record), Nx), float)  # Snapshots for the solutions

# %% First time step ---------------------------------------------------------|
for i in range(1, Nx-1):  # ensures the boundaries stay at zero
    u[i] = um[i] - 0.25*beta*(um[i+1]**2 - um[i-1]**2)  # Mind the 0.25

# %% Main loop ---------------------------------------------------------------|
ii = 0  # subplot counter
for tt, t_sec in enumerate(t):

    # Calculate the new values of u^{j+1}
    for i in range(1, Nx-1):
        up[i] = um[i] - .5*beta*(u[i+1]**2 - u[i-1]**2) + beta**2/4*(
            ( u[i]+u[i+1] )*( (u[i+1])**2-(u[i])**2 ) + (
              u[i]+u[i-1] )*( (u[i-1])**2-(u[i])**2 ) )          

    # Make plots at the given times
    if abs(t_sec-t_record[ii]) < delta:
        u_record[ii, :] = up
        if ii < len(t_record)-1:
            ii += 1
        else:
            ii = 0  # reset to zero, won't hit it again

    if tt % 10 == 0:  # we skip every 10 steps to speed up
        # continue  # comment out to see the animation
        plt.figure(1)
        plt.clf()
        # Make a plot
        plt.plot(x, up, color = 'black')   # 32 filled contours
        plt.xlabel('$x$')
        plt.ylabel('$u$')
        plt.xlim([0., Lx])
        # plt.ylim([-5e-4, 5e-4])
        plt.title('Lax-Wendroff $t = ${0:.4f}'.format(t_sec))
        plt.grid()
        if tt ==0 :
            plt.savefig('Q3_LW_t=0.png', dpi = 200)
        if tt == 50 : 
            plt.savefig('Q3_LW_t=0.5.png', dpi = 200)
        if tt == 100: 
            plt.savefig('Q3_LW_t=1.png', dpi = 200)
        if tt == 150: 
            plt.savefig('Q3_LW_t=1.5.png', dpi = 200)
        plt.draw()
        plt.pause(0.01)

    u = up
    um = u  # shift the axes in time


# %% Last plot ---------------------------------------------------------------|
plt.figure(2)
for ii in range(len(t_record)):
    plt.subplot(2, 2, ii+1)
    plt.plot(x/np.pi, u_record[ii, :])
    plt.xlabel("$x/\pi$")
    plt.ylabel("$y$")
    plt.xlim([0., Lx/np.pi])
    # plt.ylim([-5e-1, 5e-1])
    plt.title("$u$ at $t=${0:.4}".format(t_record[ii]))
    plt.grid()

plt.tight_layout()
plt.savefig('Lab08-Q3.pdf')
plt.show()
