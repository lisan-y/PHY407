"""
PHY407 Lab 10
Question 1

Generates random coordinates phi & theta.
Plots water/land distribution of Earth.

Author: Lisa Nasu-Yu, Nov 2021
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

S = 15
L = 15
T = 15

plt.rc('font', size = S)
plt.rc('axes', titlesize = T)
plt.rc('axes', labelsize = S)
plt.rc('xtick', labelsize = S)
plt.rc('ytick', labelsize = S)
plt.rc('legend', fontsize = L)
plt.rc('figure', titlesize = S)

# Define functions to generate random phi and theta in radians
def phi():
    z = np.random.random()
    return 2*np.pi*z


def theta():
    z = np.random.random()
    return np.arccos(1-2*z)


# Define function to randomly sample N points from Earth.npy and estimate land fraction
def land_distribution(N):
    longitude = [180 - 180 * phi() / np.pi for i in range(N)]
    latitude = [90 - 180 * theta() / np.pi for i in range(N)]

    # create nearest interpolator
    interp = RegularGridInterpolator((lon_closed, lat_array), data_closed, method='nearest')
    points = np.array([[longitude[i], latitude[i]] for i in range(N)])
    nearest_values = np.array(interp(points))

    # find indices of water and land components
    water_ind = np.where(np.abs(nearest_values) < 1e-15)[0]
    land_ind = np.where(np.abs(nearest_values - 1) < 1e-15)[0]
    water_lon = [longitude[int(i)] for i in water_ind]
    water_lat = [latitude[int(i)] for i in water_ind]
    land_lon = [longitude[int(i)] for i in land_ind]
    land_lat = [latitude[int(i)] for i in land_ind]

    # calculate land fraction
    land_frac = len(land_ind) / (len(land_ind) + len(water_ind))

    return land_frac, water_lon, water_lat, land_lon, land_lat


# Load Earth data
loaded = np.load('Earth.npz')
data = loaded['data']
lon_array = loaded['lon']
lat_array = loaded['lat']

# Part a
# Plot 500 random points in latitude/longitude
N = 5000
# convert to longitude/latitude.
# phi=0 corresponds to 180, theta=0 to 90
longitude = [180 - 180*phi()/np.pi for i in range(N)]
latitude = [90 - 180*theta()/np.pi for i in range(N)]
plt.figure()
plt.title(r'Random Positions')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.plot(longitude, latitude, '.')
plt.tight_layout()
plt.savefig('random.pdf')
plt.show()

# Part c
# Calculate land fraction
land = np.shape(np.where(data == 1))[1]
dphidtheta = 1/(len(lon_array) * len(lat_array))
print('Land Fraction = ', land*dphidtheta)

# Part d

# extend grid to include up to 180 longitude
lon_closed = np.concatenate((lon_array, [180]))
# set longitude=180 data equal to longitude=-180 data
data_closed = np.concatenate((data, [data[0]]), axis=0)
lat_grid, lon_grid = np.meshgrid(lat_array, lon_closed)

# calculate land frac for various N
N = [50, 500, 5000]
land_fracs = []
for n in N:
    land_frac, _,_,_,_ = land_distribution(n)
    land_fracs.append(land_frac)

# calc land fraction and plot distribution for N=50000
N = 50000
land_frac, water_lon, water_lat, land_lon, land_lat = land_distribution(N)
land_fracs.append(land_frac)

print('Land Fractions for N=50,500,5000,50000: {}'.format(land_fracs))

plt.figure()
plt.title('Land/Water Distribution on Earth')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.xlim(-180,180)
plt.ylim(-90,90)
plt.plot(water_lon, water_lat, '.', label='Water')
plt.plot(land_lon, land_lat, '.', color='green', label='Land')
plt.legend()
plt.tight_layout()
# plt.savefig('earth.pdf')
plt.show()

