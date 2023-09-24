#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 11:12:49 2021

@author: emmajarvis
"""

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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
    

#Functions from Lab Manual:

def get_tau_step():
    """ 
    Calculate how far a photon travels before it gets scattered.
    OUT: 
        optical depth traveled 
    """
    delta_tau = -np.log(np.random.random())
    return delta_tau

def emit_photon(tau_max):
    """ 
    Emit a photon from the stellar core.
    IN: 
        tau max is max optical depth
    OUT:
        tau: optical depth at which the photon is created
        mu: directional cosine of the photon emitted
    """
    tau = tau_max
    delta_tau = get_tau_step()
    mu = np.random.random()
    return tau - delta_tau*mu, mu

def scatter_photon(tau):
    """ 
    Scatter a photon.
    IN: 
        tau, optical depth of the atmosphere
    OUT:
        tau: new optical depth
        mu: directional cosine of the photon scattered
    """
    delta_tau = get_tau_step()
    mu = 2*np.random.random()-1 # sample mu uniformly from -1 to 1
    tau = tau - delta_tau*mu
    return tau, mu

def photon_path(tau_max):
    
    tau, mu = emit_photon(tau_max)
   
    N = 0 #number of times it scatted
    while tau > 0:
        tau, mu = scatter_photon(tau)
        N += 1
        if tau > tau_max:
            tau, mu = emit_photon(tau_max)
            N = 0
    
    return mu, N

#%% Obtain the mus and N for 100 000 photons
all_mus = []
all_Ns = []



for j in tqdm(range(100000)):
    mu, N = photon_path(10)
    all_mus.append(mu)
    all_Ns.append(N)
  
     
all_mus = np.array(all_mus)
all_Ns = np.array(all_Ns)
#%%


Ns = all_Ns

hist, binedges = np.histogram(all_mus, bins = 21)
hist = hist/max(hist)

#obtain the midpoints of the bins
mu_mid = []


for i in range(len(binedges)-1):
    
    mu_av = (binedges[i+1]+binedges[i])/2
     
    mu_mid.append(mu_av)

mu_mid = np.array(mu_mid)

# print(mu_mid, N_mid)

plt.figure()
plt.bar(binedges[:len(binedges)-1], hist, 
        width = 0.95*np.abs(binedges[1]-binedges[0]), 
        alpha = 0.5, color = 'purple')
griddy()
plt.xlabel(r'$\mu$')
plt.ylabel(r'$N/N_1$')
plt.title(r'$\tau = 10$')
# plt.savefig('/Users/emmajarvis/Desktop/Q2b_histogram_mu.png', dpi = 200)
plt.show()


def intensity(mu, a, b):
    return a*mu + b

#Intensity
I = hist/mu_mid

popt, pcov = curve_fit(intensity, mu_mid, I)

print('curve-fitted coefficients: ', popt)
print('curve-fitted coefficients, error: ', pcov)
print('expected coefficient: ', [0.6, 0.4])

I_analytic = intensity(all_mus, 0.6, 0.4)
I_fit = intensity(mu_mid, popt[0], popt[1])


plt.figure()
plt.plot(all_mus, I_analytic, label = 'Analytic Solution', 
         color = 'dodgerblue', linestyle = '--', linewidth = 1)
plt.plot(mu_mid, I_fit, label = 'Fitted solution', color = 'purple')
griddy()
plt.xlabel(r'$\mu$')
plt.ylabel(r'$I(\mu)/I_1$')
plt.title(r'$\tau = 10$')
plt.legend()
# plt.savefig('/Users/emmajarvis/Desktop/Q2b_intensity.png', dpi = 200)
plt.show()

# print('part b complete')
#%% Repeat to get mus and N but with tau = 1e-4

all_mus_c = []
all_Ns_c = []



for j in tqdm(range(100000)):
    mu, N = photon_path(1e-4)
    all_mus_c.append(mu)
    all_Ns_c.append(N)
  
     
all_mus_c = np.array(all_mus_c)
all_Ns_c = np.array(all_Ns_c)
#%% Repeat the histogram and curve fitting for the new tau
Ns_c = all_Ns_c

hist_c, binedges_c = np.histogram(all_mus_c, bins = 21)
hist_c = hist_c/max(hist_c)

#obtain the midpoints of the bins
mu_mid_c = []


for i in range(len(binedges_c)-1):
    
    mu_av_c = (binedges_c[i+1]+binedges_c[i])/2
     
    mu_mid_c.append(mu_av_c)

mu_mid_c = np.array(mu_mid_c)

# print(mu_mid, N_mid)

plt.figure()
plt.bar(binedges_c[:len(binedges_c)-1], hist_c, 
        width = 0.95*np.abs(binedges_c[1]-binedges_c[0]), 
        alpha = 0.5, color = 'green')
# plt.bar(all_mus, all_Ns)
griddy()
plt.xlabel(r'$\mu$')
plt.ylabel(r'$N/N_1$')
plt.title(r'$\tau = 10^{-4}$')
# plt.savefig('/Users/emmajarvis/Desktop/Q2c_histogram_mu.png', dpi = 200)
plt.show()




def intensity(mu, a, b):
    return (a*mu + b)

#Intensity
I_c = hist_c/mu_mid_c

popt_c, pcov_c = curve_fit(intensity, mu_mid_c, I_c)

print('curve-fitted coefficients: ', popt_c)
print('expected coefficient: ', [0.6, 0.4])

I_analytic_c = intensity(all_mus_c, 0.6, 0.4)
I_fit_c = intensity(mu_mid_c, popt_c[0], popt_c[1])



plt.figure()
plt.plot(all_mus_c, I_analytic_c, label = 'Analytic Solution', 
         color = 'dodgerblue', linestyle = '--', linewidth = 1)
plt.plot(mu_mid_c, I_fit_c, label = 'Fitted solution', color = 'green')
griddy()

plt.xlabel(r'$\mu$')
plt.ylabel(r'$I(\mu)/I_1$')
plt.title(r'$\tau = 10^{-4}$')
plt.legend()
# plt.savefig('/Users/emmajarvis/Desktop/Q2c_intensity.png', dpi = 200)
plt.show()


