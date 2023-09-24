#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 12:56:07 2021

@author: emmajarvis
"""

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

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
#%% Q3a    
   
N = 10000
a = 0
b = 1

def f(x):
    return x**-0.5/(np.exp(x)+1)

#%% Mean Value Method integral values

all_I_mvm = []
for j in tqdm(range(100)):
    k = 0
    for i in range(N):
        x = (b-a)*np.random.random()
        k += f(x)
        I = k * (b-a) / N
    all_I_mvm.append(I)
    
#%%  Plot Histogram of integral using mean value method 

plt.figure()
plt.hist(all_I_mvm, 10, range=[0.8,0.88], color = 'black', alpha = 0.5)
plt.xlabel('I')
plt.ylabel('N')
plt.title(r'Mean Value Method, $I = \int^1_0 \frac{x^{-1/2}}{1+e^x}dx$')
griddy()
plt.tight_layout()
# plt.savefig('/Users/emmajarvis/Desktop/Q3a_mean_value.png', dpi = 200)
plt.show()
    
#%% Importance Sampling integral values

all_I_is = []
for j in tqdm(range(100)):
    k = 0
    for i in range(N):
        z = (b-a)*np.random.random()
        x = z**0.5
        k += np.exp(x)*f(x)
        I = k/N
    all_I_is.append(I)

#%%  Plot Histogram of integral using importance sampling method 

plt.figure()
plt.hist(all_I_is, 10, range=[0.8,0.88], color = 'red', alpha = 0.5)
plt.xlabel('I')
plt.ylabel('N')
plt.title(r'Importance Sampling Method, $I = \int^1_0 \frac{x^{-1/2}}{1+e^x}dx$')
griddy()
plt.tight_layout()
# plt.savefig('/Users/emmajarvis/Desktop/Q3a_importance_sampling.png', dpi = 200)
plt.show()

#%% Q3b

a = 0
b = 10
N = 10000


def f(x):
    return np.exp(-2*np.abs(x-5))
    
#%% Mean Value Method integral values

all_I_mvm = []
for j in tqdm(range(100)):
    k = 0
    for i in range(N):
        x = (b-a)*np.random.random()
        k += f(x)
        I = k * (b-a) / N
    all_I_mvm.append(I)

#%%  Plot Histogram of integral using mean value method 

plt.figure()
plt.hist(all_I_mvm, 10, range=[0.95,1.05], color = 'black', alpha = 0.5)
plt.xlabel('I')
plt.ylabel('N')
plt.title(r'Mean Value Method, $I = \int^{10}_0 e^{-2|x-5|}dx$')
griddy()
plt.tight_layout()
# plt.savefig('/Users/emmajarvis/Desktop/Q3b_mean_value.png', dpi = 200)
plt.show()

#%% Importance Sampling integral values

all_I_is = []
for j in tqdm(range(100)):
    k = 0
    for i in range(N):
        x = np.random.normal(5, 1)
        k +=  3 * f(x)
        I = k/N
    all_I_is.append(I)
    
#%%  Plot Histogram of integral using importance sampling method 

plt.figure()
plt.hist(all_I_is, 10, range=[0.95,1.05], color = 'red', alpha = 0.5)
plt.xlabel('I')
plt.ylabel('N')
plt.title(r'Importance Sampling Method, $I = \int^{10}_0 e^{-2|x-5|}dx$')
griddy()
plt.tight_layout()
# plt.savefig('/Users/emmajarvis/Desktop/Q3b_importance_sampling.png', dpi = 200)
plt.show()




