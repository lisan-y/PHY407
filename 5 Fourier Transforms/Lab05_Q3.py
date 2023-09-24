#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 12:41:18 2021

@author: Emma Jarvis

Lab05 Q3 - all parts.
"""

import numpy as np
import matplotlib.pyplot as plt

# Plotting parameters
plt.rcParams['figure.figsize'] = (10,6)
plt.rcParams['font.family']='serif'

plt.rcParams['mathtext.fontset'] = 'cm'

S = 20
L = 15
T = 20

plt.rc('font', size = S)
plt.rc('axes', titlesize = T)
plt.rc('axes', labelsize = S)
plt.rc('xtick', labelsize = S)
plt.rc('ytick', labelsize = S)
plt.rc('legend', fontsize = L)
plt.rc('figure', titlesize = S)

def griddy():
    plt.minorticks_on()
    plt.grid(color='black',
             which='major',
             linestyle=":",
             linewidth='0.4',
             )
    plt.grid(color='black',
             which='minor',
             linestyle=":",
             linewidth='0.2',
             )


directory = ''  # directory containing txt files

#Load files
SLP = np.loadtxt(directory + 'SLP.txt')
Longitude = np.loadtxt(directory + 'lon.txt')
Times = np.loadtxt(directory + 'times.txt')



"""
Q3a
"""
cmap = 'twilight'

#Make contour plot of original SLP data
plt.contourf(Longitude, Times, SLP, cmap = cmap)
plt.xlabel('longitude (degrees)')
plt.ylabel('days since Jan. 1 2015')
plt.title('SLP anomaly (hPa)')
plt.colorbar(label = 'Sea Level Pressure (hPa)')
# plt.savefig('/Users/emmajarvis/Desktop/SLP_data_no_fft.png', dpi = 300)
plt.show()



rfft_SLP = np.fft.rfft(SLP) #take fft of SLP data

rfft_3 = rfft_SLP[:,3] #take slice corresponding to wavenumber 3
wavenumber_3 = np.zeros(np.shape(rfft_SLP), dtype = complex) 
wavenumber_3[:,3] = rfft_3
inverse_3 = np.fft.irfft(wavenumber_3) #take inverse, get back into position space


rfft_5 = rfft_SLP[:,5] #take slice corresponding to wavenumber 5
wavenumber_5 = np.zeros(np.shape(rfft_SLP), dtype = complex) 
wavenumber_5[:,5] = rfft_5
inverse_5 = np.fft.irfft(wavenumber_5) #take inverse, get back into position space

#Contour plot of m = 3 wavenumber
plt.contourf(inverse_3, cmap = cmap)
plt.colorbar(label = 'Sea Level Pressure (hPa)')
plt.xlabel('longitude (degrees)')
plt.ylabel('days since Jan. 1 2015')
plt.title('SLP anomaly (hPa), m = 3')
# plt.savefig('/Users/emmajarvis/Desktop/SLP_data_m=3.png', dpi = 300)
plt.show()

#Contour plot of m = 5 wavenumber
plt.contourf(inverse_5, cmap = cmap)
plt.colorbar(label = 'Sea Level Pressure (hPa)')
plt.xlabel('longitude (degrees)')
plt.ylabel('days since Jan. 1 2015')
plt.title('SLP anomaly (hPa), m = 5')
# plt.savefig('/Users/emmajarvis/Desktop/SLP_data_m=5.png', dpi = 300)
plt.show()

"""
Q3b
"""

lat = 50
r_E = 6365000
circ = 2*np.pi*r_E*np.cos(lat)
print(' \n The circumference of the latitude is ', circ, 'm')





    
