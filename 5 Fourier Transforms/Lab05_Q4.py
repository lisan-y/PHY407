#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 10:22:56 2021

@author: emmajarvis
"""
import numpy as np
import matplotlib.pyplot as plt

"""
Q4a 
"""

#Load the data
# blur = np.loadtxt('/Users/emmajarvis/Desktop/cpresources/Blur.txt')

blur = np.loadtxt('Blur.txt')

#Plot the original blurry image
plt.imshow(blur, cmap = 'gray')
plt.tight_layout()
# plt.savefig('/Users/emmajarvis/Desktop/blurry_image.png', dpi = 300)
plt.show()

"""
Q4b 
"""


def gaussian(x, y, sigma):
    """
    Obtain the gaussian point spread function at a point (x, y) with a 
    width of sigma.

    Parameters
    ----------
    x : int
        x value of a point.
    y : int
        y value of a point.
    sigma : int
        width of point spread function.

    Returns
    -------
    float
        point spread function at point (x,y).

    """
    return np.exp(-(x**2+y**2)/(2*sigma**2)) #equation defined in Newman 7.9


#Value of width, given in equation
sigma = 25

def get_gauss_image():
    """

    Returns
    -------
    gaussian_image : 2D array
        2D array corresponding to the value obtained using gaussian function
        for each point in a (1024,1024) array.

    """
    gaussian_image = np.zeros(np.shape(blur), dtype = complex) #initialize array of zeros to hold the point spread function
    for i in range(-512, 512, 1):
        for j in range(-512, 512, 1):
            gaussian_image[i][j] = gaussian(i, j, sigma) #apply gaussian function to each point
    return gaussian_image
 
#obtain point spread function using get_gauss_image and plot it using imshow    
gaussian_image = get_gauss_image()
plt.imshow(abs(gaussian_image), cmap = 'gray')
plt.tight_layout()
# plt.savefig('/Users/emmajarvis/Desktop/gaussian_image.png', dpi = 300)
plt.show()


"""
Q4c 
"""

def unblur():
    """
    Unblur an image.

    Returns
    -------
    Plot of unblurry image.

    """
    
    image_fft = np.fft.rfft2(blur) #take fourier transform of the blurry image
    
    gaussian_fft = np.fft.rfft2(gaussian_image) #take the fourier transform of the point spread function
    
    unblurred = np.zeros(image_fft.shape, dtype = complex) # initializing an array of zeros to become the unblurred photo
    
    #loop over each point
    for i in range(gaussian_fft.shape[0]):
        for j in range(gaussian_fft.shape[1]):
            if gaussian_fft[i][j] > 1e-5: #avoid divide by zero error 
                unblurred[i][j] = image_fft[i][j]/gaussian_fft[i][j] #divide image by point spread function
            else:
                unblurred[i][j] = image_fft[i][j] #avoid divide by zero error by not dividing by the gaussian
                
                
    unblurred_image = np.fft.irfft2(unblurred) #inverse Fourier transform
    plt.imshow(abs(unblurred_image), cmap ='gray') #displat image

#Run the unblur funtion:    
unblur()
# plt.savefig('/Users/emmajarvis/Desktop/unblurred_image.png', dpi = 300)
plt.show()



