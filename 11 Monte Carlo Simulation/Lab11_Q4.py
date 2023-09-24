#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 12:35:26 2021

Author: Emma Jarvis

Lab 11 Q4: Ising Model

Code adapted from L11-Ising1D-start.py
"""


# import modules
import numpy as np
from random import random, randrange
import matplotlib.pyplot as plt
import matplotlib as mpl

# %matplotlib qt

#Plotting Parameters:
plt.rcParams['figure.figsize'] = (10,8)
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
    
    
#%% parts a, b

def energyfunction(J_, dipoles):
    """ function to calculate energy """
    adjacentx = np.zeros([1,20])
    adjacenty = np.zeros(20)
    
    for i in range(19):
          adjacentx += dipoles[:, i]*dipoles[:,i+1]
          adjacenty += dipoles[i,:]*dipoles[i+1, :]
    total = J_*(np.sum(adjacentx[0])+np.sum(adjacenty))
    return total


def acceptance(Ej, Ei, kT):
    """ Function for acceptance probability """
    result = False
    
    if Ei-Ej <0:
        return True
    if random() < np.exp(Ej-Ei)/kT:
        result = True
    
    return result  # result is True of False


def ising(T, animation = False, N = int(1e5)):
    """ 
    use Markov chain Monte Carlo Simulation for Ising model
    
    T = temperature
    N = number of steps
    
    """
    J = 1.0
    kB = 1.0
    num_dipoles = 20

    # generate array of dipoles and initialize diagnostic quantities
    dipoles = np.ones([num_dipoles, num_dipoles], int)  
    for n in range(num_dipoles):
        for m in range(num_dipoles):
            flip = random()
            if flip < 0.5:
                dipoles[n][m] *= -1
    
    #lists to store energy and magnetization values        
    energy = []  
    magnet = []

    #loop over all steps
    for i in range(N):
        
        if i%(N/50) == 0:
            print('Step '+str(i)+'/'+str(N))
        if i%(N/2) == 0:
            print('halfway there...')    
        
        #randomly pick one of the dipoles
        picked1 = randrange(num_dipoles)
        picked2 = randrange(num_dipoles) 
        
        #old energy
        E = energyfunction(J, dipoles)
        energy.append(E)
        magnet.append(np.sum(dipoles))
        
        # propose to flip the victim
        dipoles[picked1][picked2] *= -1  
        
        # compute Energy of proposed new state
        Enew = energyfunction(J, dipoles)  
        Mnew = np.sum(dipoles)
        
        # calculate acceptance probability 
        accepted = acceptance(Enew, E, kB*T)
    
        # store energy and magnetization
        if accepted:
            energy.append(Enew)
            magnet.append(Mnew)
            
        else:
            energy.append(E)
            magnet.append(np.sum(dipoles))
            dipoles[picked1][picked2] *= -1
        
        #make the animation
        if animation == True:
            if i%(N/50) ==0:
                plt.clf()
                
                cmap = mpl.cm.coolwarm
                bounds = [-1, 1]
                norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')
                plt.imshow(dipoles, vmin = -1, vmax = 1, cmap = cmap)
                plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), 
                             label = 'Spin', extend = 'both', 
                             spacing = 'proportional',
                             extendfrac ='auto',
                             aspect = 5)
                plt.title('Step Number = '+str(i)+ ', Temperature = '+str(T))
                plt.draw()
                plt.pause(0.01)
                # if i == 30000 or i == 60000 or i == 90000:
                    # plt.savefig('/Users/emmajarvis/Desktop/Q4e_T='+str(T)+'step='+str(i)+'.png', dpi = 200)
    
    return energy, magnet, dipoles
#%% part c, d

#plot energy and magnetization for 1 million steps
energy, magnet, dipoles = ising(1.0, N = int(1e6))

plt.figure()
plt.plot(energy, linestyle = '-', marker = '.', color = 'black', ms = 2)
plt.ylabel('Energy')
plt.xlabel('Number of Monte Carlo Steps')
griddy()
# plt.savefig('/Users/emmajarvis/Desktop/Q4_c_energy.png')
plt.show()


plt.figure()
plt.plot(magnet, linestyle = '-', marker = '.', color = 'black', ms = 2)
plt.ylabel('Megnetization')
plt.xlabel('Number of Monte Carlo Steps')
griddy()
# plt.savefig('/Users/emmajarvis/Desktop/Q4_c_magnetization.png')
plt.show()

#%% part e

#make the animation for temperatures of 1, 2 and 3
energy1, magnet1, dipoles1 = ising(1.0, animation = True)
energy2, magnet2, dipoles2 = ising(2.0, animation = True)
energy3, magnet3, dipoles2 = ising(3.0, animation = True)
