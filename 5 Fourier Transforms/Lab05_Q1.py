"""
PHY407 Lab 05
Question 1

a) Exercise 7.1, Newman
Plots coefficients of the discrete Fourier transform of several periodic functions.

b) Exercise 7.8, Newman
Computes and plots diffraction grating pattern by taking a Fourier Transform

Author: Lisa Nasu-Yu, Oct 2021
"""

import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import matplotlib as mpl


# plot settings
S = 15
L = 15
T = 15

plt.rc('font', size=S)
plt.rc('axes', titlesize=T)
plt.rc('axes', labelsize=S)
plt.rc('xtick', labelsize=S)
plt.rc('ytick', labelsize=S)
plt.rc('legend', fontsize=L)
plt.rc('figure', titlesize=S)

# set plot fonts. runs slower

# mpl.rcParams['legend.frameon'] = False
# mpl.rcParams['figure.autolayout'] = True
#
# plt.rcParams.update({"text.usetex":True, "font.family": 'sans-serif', 'font.sans-serif': ['Helvetica']})
# plt.rcParams.update({"text.usetex":True, "font.family": 'serif', 'font.serif': ['Palatino']})


# Q1 a), 7.1 a)


def dft(y):
    """
    Calculates a DFT the slow way.

    Written by Mark Newman <mejn@umich.edu>, June 4, 2011
    :param y: [array] sample points
    :return: array of values from fourier transform
    """
    N = len(y)
    c = np.zeros(N//2+1, complex)
    for k in range(N//2+1):
        for n in range(N):
            c[k] += y[n] * cm.exp(-2j * np.pi * k * n / N)
    return c


def plot(y, c, title):
    """
    Plots a set of sample points and |c_k| from its fourier transform
    :param y: [array] sample points
    :param c: [array] fourier transformed values
    :param title: [str] title of the function plot
    :return: 2 plots
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    # plot function
    ax1.set_title(title)
    ax1.set_ylabel(r'y_n')
    ax1.set_xlabel(r'n')
    ax1.plot(y)

    # plot fourier transform
    ax2.set_title('Fourier Transform')
    ax2.set_xlabel(r'k')
    ax2.set_ylabel(r'$|c_k|$')
    ax2.bar(range(len(c)), np.abs(c), width=1.5)

    # show cropped plot of fourier transform for better viewing
    ins = ax2.inset_axes([0.4, 0.4, 0.6, 0.6])
    ins.set_xlabel(r'k')
    ins.set_ylabel(r'$|c_k|$')
    ins.set_xlim(-1, 25)
    ins.bar(range(len(c)), np.abs(c), width=1.1)
    plt.tight_layout()
    plt.savefig(title+'.pdf')
    plt.show()


# Set number of samples
N = 1000
# Define sample points
y_a = np.zeros(N, complex)
y_a[:N // 2] = 1.
# fourier transform
c_a = dft(y_a)
# # plot
plot(y_a, c_a, 'Square Wave')

# Q1 a), 7.1 b)
# Define sample points
y_b = np.linspace(0, N, 1000, complex)
# fourier transform
c_b = dft(y_b)
# plot
plot(y_b, c_b, 'Sawtooth Wave')

# Q1 a), 7.1 c)
# Define sample points
y_c = [np.sin(np.pi * n / N) * np.sin(20 * np.pi * n / N) for n in range(N)]
# fourier transform
c_c = dft(y_c)
# plot
plot(y_c, c_c, 'Modulated Sine Wave')


# Q1 b), 7.8

def y(u):
    """
    Computes y_n, with padding points=0 at end of array.
    :param u: [float] sampling point
    :return: [float] value at sampling point
    """
    # set padding sample points at end of array
    if np.abs(u) > w / 2:
        return 0
    # compute y
    else:
        return np.sqrt(np.sin((np.pi / w_original) * u) ** 2)


# Define constants
w_original = 20e-6  # slit width [m]
w = w_original * 10  # grating width [m]
W = 10 * w  # [m]

wavelength = 500e-9  # [m]
f = 1  # focal length [nm]
screen = 0.1  # width of screen [m]
dx = wavelength * f / w  # [m]
N = 1000  # number of sampling points

# define sample points
u = np.arange(0, N) * W / N - W / 2
# compute y_n from sample points
y_n = np.array([y(i) for i in u])
# apply fourier transform
y_ft = np.fft.rfft(y_n)
# compute FT
I = (W ** 2 / N ** 2) * np.abs(y_ft) ** 2

# define x range
x_points = np.arange(-0.5, 0.5, dx)
# Set second half of full length N array of I as complex conj of reflected first half
# Note that we flip the first and second half of I because we want the complex conjugate of the FT.
I_full = np.concatenate((I[len(x_points) // 2:0:-1], I[:len(x_points) // 2 + 1]))

# plot
plt.figure()
plt.title('Diffraction pattern')
plt.ylabel(r'Intensity [$W/m^2$]')
plt.xlabel('x [m]')
plt.plot(x_points, I_full)
plt.tight_layout()
plt.savefig('diffraction.pdf')
plt.show()


