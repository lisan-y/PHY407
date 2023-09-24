# %load laplace.py
# Newman's laplace.py
from numpy import empty, zeros, max, min, log10
from pylab import imshow, gray, show
import matplotlib.pyplot as plt

# Constants
M = 100         # Grid squares on a side
V = 1.0         # Voltage at top wall
target = 1e-6   # Target accuracy
w=0.9

# Create arrays to hold potential values
phi = zeros([M+1, M+1], float)
phi[0, :] = V
phiprime = empty([M+1, M+1], float)

# Main loop
delta = 1.0
while delta > target:

    # Calculate new values of the potential
    for i in range(M):
        for j in range(M):
            #if i == 0 or i == M or j == 0 or j == M:
            if 1==0:
                phiprime[i, j] = phi[i, j]
            else:
                phiprime[i, j] = (1+w)*(phi[i+1, j] + phi[i-1, j]
                                  + phi[i, j+1] + phi[i, j-1])/4 - w*phi[i,j]

    # Make a plot
    plt.clf()
    plt.imshow(phi)
    plt.gray()
    plt.colorbar()
    plt.draw()
    plt.pause(0.1)

    # Calculate maximum difference from old values
    delta = max(abs(phi-phiprime))
    print(delta)

    # Swap the two arrays around
    phi, phiprime = phiprime, phi

# Make a plot
plt.imshow(phi)
# plt.scale('log')
plt.gray()
plt.colorbar()
plt.show()