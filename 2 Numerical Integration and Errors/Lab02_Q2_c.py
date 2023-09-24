# Qc

S = 20
L = 20
T = 20

plt.rc('font', size = S)
plt.rc('axes', titlesize = T)
plt.rc('axes', labelsize = S)
plt.rc('xtick', labelsize = S)
plt.rc('ytick', labelsize = S)
plt.rc('legend', fontsize = L)
plt.rc('figure', titlesize = S)


def u(n, r, theta):
    # set R (boundary radius) and c (wave phase speed) as 1
    # set z_mn = z_3,2
    R = 1
    c = 1
    z_mn = 11.620
    # note, the second cosine term, cos(c*z_m*t/a) = 1 for t = 0
    return bessel(n, z_mn * r / R) * np.cos(n * theta)



'''
relies on a value that may not exist at every point
taylor expansion, what to do if lower order terms are zero?
based on third order term that vanishes at boundary
so have to consider next order term
explain why we don't have the third term. 
why might it not work

certain values vanish in our equation

solution: shift all data by mean
'''

from matplotlib import cm
from matplotlib.ticker import LinearLocator


# u[R > R_0] = 0.0

# 3 angles

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
dtheta = 0.01
X = np.cos(np.arange(0, 2 * np.pi, dtheta))
Y = np.sin(np.arange(0, 2 * np.pi, dtheta))
X, Y = np.meshgrid(X, Y)
r = np.sqrt(X ** 2 + Y ** 2)
Z = u(2, r, np.arange(0, 2 * np.pi, dtheta))
Z[r > 1] = 0.0

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# set ticks
ax.set_xticks([-1, 0, 1])
ax.set_yticks([-1, 0, 1])
ax.set_zticks([-0.5, 0, 0.5])

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
