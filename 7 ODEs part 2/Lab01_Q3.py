#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 19:50:47 2021

@author: emmajarvis
"""

"""
Code that histograms N random samples into M linearly spaced bins.

"""

import numpy as np
from time import time
import matplotlib.pyplot as plt


"""
Plotting Information
"""
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
plt.rcParams['figure.figsize'] = (10,6)
plt.rcParams['font.family']='serif'

plt.rcParams['mathtext.fontset'] = 'cm'

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


#3a
# Pseudocode
# Import numpy and matplotlib.pyplot and time
# set the minimum bin
# set the maximum bin
# generate an array of M linearly spaced bins between -5 and 5
# generate N samples using numpy.random.randn
# sort N samples into the M bins


"""
FUNCTION NAME:
make_hist_no_np

PURPOSE:
To generate a histogram from scratch given an array of samples and number of bins.

PARAMETERS:
samples: An array of N random samples.
M: number of linearly spaced bins.
plot: if plot == True this function generates a plot of the histogram. Default is false.

RETURNS:
None
if plot == True, a histogram is generated.
"""
def make_hist_no_np(samples, M, plot = False):
    hmin = -5 # lowest bin in histogram to plot
    hmax = 5 #highest bin in histogram to plot
    bins = np.linspace(hmin, hmax, M) #Create bins
    spacing = bins[1]-bins[0]
    hist = np.array([np.where(np.logical_and(samples>i, samples<i+spacing))[0].size for i in bins]) #y-axis of histogram
    
    if plot == True:
        hist = hist/len(samples) #Normalize
        plt.figure()
        plt.bar(bins, hist, color = 'deepskyblue', alpha = 0.5)
        griddy()
        plt.title('Histogram without numpy.histogram, N = '+str(len(samples)))
        plt.ylabel('Probability')
        plt.xlabel('Sample Number')
        plt.tight_layout()
        plt.savefig('/Users/emmajarvis/Desktop/Histogram_no_numpy_N_'+str(len(samples))+'.png')
        plt.show()
    

#3b
"""
FUNCTION NAME:
make_hist_with_np

PURPOSE:
To generate a histogram using np.histogram given an array of samples and number of bins.

PARAMETERS:
samples: An array of N random samples.
M: number of linearly spaced bins.
plot: if plot == True this function generates a plot of the histogram. Default is false.

RETURNS:
None
if plot == True, a histogram is generated.
"""
def make_hist_with_np(samples, M, plot = False):
    hist, binedges = np.histogram(samples, bins = np.linspace(-5, 5, M))
    hist = hist/len(samples) #Normalize
    if plot == True:
        plt.figure()
        griddy()
        plt.bar(binedges[:len(binedges)-1], hist, color = 'mediumpurple', alpha = 0.5)
        
        plt.ylabel('Probability')
        plt.xlabel('Sample Number')
        plt.title('Histogram with numpy.histogram, N = '+str(len(samples)))
        plt.tight_layout()
        plt.savefig('/Users/emmajarvis/Desktop/Histogram_with_numpy_N_'+str(len(samples))+'.png')
        plt.show()
        
"""
FUNCTION NAME:
plot_N_t

PURPOSE:
To time how long it takes make_hist_no_no and make_hist_with_np to generate data for a histogram
(not including the tine it takes to plot) and then plot these time values as a function of the
number of samples.

PARAMETERS:
N: An array the number of samples to be generated
M: number of linearly spaced bins.
plot: if plot == True this function generates a plot of the histogram. Default is false.

RETURNS:
None
if plot == True, histograms are generated for each method and each value of N.
plots the time as a function of the number of samples, N.
"""
def plot_N_t(N, M, plot = False):
    time_no_np = [] #list to include time values for the histogram generated without np.histogram
    time_with_np = [] #list to include time values for the histogram generated using np.histogram
   
    for num in N:
        samples = np.random.randn(num) #generates the array of num random samples for every num in N.
        
        start = time()
        make_hist_no_np(samples, M, plot) #timing the histogram generated without np.histogram
        end = time()
        diff = end-start
        time_no_np.append(diff) #add this time value to the list
        
        start2 = time()
        make_hist_with_np(samples, M, plot) #timing the histogram generated using np.histogram
        end2 = time()
        diff2 = end2-start2
        time_with_np.append(diff2) #add this time value to the list
    
    print('Times without numpy = ', time_no_np) #prints the lists of time values  
    print('Times with numpy = ', time_with_np)

    #PLOT
    plt.figure()
    plt.plot(N, time_no_np, label = 'without using numpy.histogram', color = 'deepskyblue', linewidth = 3)
    plt.semilogy(N, time_with_np, label = 'using numpy.histogram', color = 'mediumpurple', linewidth = 3)
    plt.legend()
    plt.ylabel('Time (s)')
    plt.xlabel('Number of Samples')
    plt.title('Time Taken for Histogram Function to Run')
    griddy()
    plt.tight_layout()
    plt.savefig('Samples_vs_time.png')
    plt.show()

"""
Perform the function that times the histograms and generates a plot of N vs. t.
"""
def main():
    N = [10, 100, 1000, 10000, 100000, 1000000] # Number of samples
    M = 1000 # Number of bins    
    plot_N_t(N, M)
if __name__ == '__main__':
    main()
