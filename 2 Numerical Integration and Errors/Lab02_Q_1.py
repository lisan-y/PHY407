#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 09:01:27 2021

@author: emmajarvis
"""
"""
IMPORTS
"""
import numpy as np

"""
Pseudocode that calculates the standard deviation using equations (1) and (2)
from the second lab assignment. Then, the relative errors for this two methods
are computed. The correct value is taken to be the value of the standard 
deviation computed using numpy.std. The relative error for each method is 
computed by subtracting the correct value from the other value, and then 
dividing by the correct value.
"""

# Pseudocode
# 1. Import numpy

# 2. Compute standard deviation using eq (1)
#   - Calculate the mean of the values
#   - Iterating over each value:
    #   - Subtract the mean
    #   - Square this value
    #   - Compute the sum of these values
#   - Divide the sum by 1 less than the total number of values
#   - Take the square root of this sum

# 3. Compute the standard deviation using eq (2)
#   - Sum over the square of all the values
#   - Subtract the mean*number of values from this sum
#   - Divide this value by 1 less than the total number of values
#   - Take the square root of the absolute value of this sum


# 4. Compute standard deviation using numpy.std

# 5. Compute the relative error of the values, x, to the true values, y, 
#       computed using numpy.std as (x-y)/y


"""
FUNCTION: calc_std_1

Purpose:
    Compute the standard deviation using the first method.
    
Parameter:
    values : array of values
    
Putpuse:
    std: standard deviation
"""
def calc_std_1(values):
    n = len(values)  # Total number of values
    Sum = np.sum(values)  # Sum of all the values
    Mean = Sum/n
    
    s = 0 #sum of sqares of the mean subtracted from each value
    for i in range(n):
        s += (values[i] - Mean)**2 # add the square of the value minus the mean squared
    
    s = (1/(n-1)) * s #divide by one less than the total number of vales
    std = np.sqrt(s) #compute the standard deviation by taking the square root
    
    return std
     
"""
FUNCTION: calc_std_2

Purpose:
    Compute the standard deviation using the second method.
    
Parameter:
    values : array of values
    
Putpuse:
    std: standard deviation
"""
def calc_std_2(values):
    n = len(values)
    s = np.sum((values)**2)
    
    inside_sqrt = (1/(n-1)) * (s - n * (np.mean(values))**2)
    
    if inside_sqrt < 0:
        print ('Error: negative value inside the square root')
    
    std = np.sqrt(abs(inside_sqrt))
    
    return std
   

"""
FUNCTION: relative_error

Purpose: compute the relative error between two values

Parameters: 
    values: array of values
    std_method: the standard deviation computation method to be used
    
Return: 
    rel_error: the relative error
"""     
def relative_error(values, std_method):
    correct_std = np.std(values, ddof = 1)
    other_std = std_method(values)
    
    rel_error = (other_std-correct_std)/correct_std
    
    return rel_error


#Q1.b

print('\n', 'Question 1.b')
values = np.loadtxt('cdata.txt')  # cdata.txt must be in same folder
std1 = calc_std_1(values)
print('The standard deviation using eq. (1) is ', std1)

rel_error_1 = relative_error(values, calc_std_1)
print('The relative error using eq. (1) is ', rel_error_1)

std2 = calc_std_2(values)
print('The standard deviation using eq. (2) is ', std2)

std_np = np.std(values, ddof = 1)
print('The standard deviation using np.std is ', std_np)

rel_error_2 = relative_error(values, calc_std_2)
print('The relative error using eq. (2) is ', rel_error_2)


#Q1.c

print('\n', 'Question 1.c')
mean1, sigma1, n1 = (0., 1., 2000)
vals1 = np.random.normal(mean1, sigma1, n1)

rel_error_1_0 = relative_error(vals1, calc_std_1)
print('The relative error using eq. (1) for mean = 0 is ', rel_error_1_0)

rel_error_2_0 = relative_error(vals1, calc_std_2)
print('The relative error using eq. (2) for mean = 0 is ', rel_error_2_0)

mean2, sigma2, n2 = (1.e7, 1., 2000)
vals2 = np.random.normal(mean2, sigma2, n2)

rel_error_1_1e7 = relative_error(vals2, calc_std_1)
print('The relative error using eq. (1) for mean = 1e7 is ', rel_error_1_1e7)

rel_error_2_1e7 = relative_error(vals2, calc_std_2)
print('The relative error using eq. (2) for mean = 1e7 is ', rel_error_2_1e7)



#Q1.d
    
"""
FUNCTION: calc_std_3

Purpose:
    Compute the standard deviation by modifying the second method.
    
Parameter:
    values : array of values
    
Putpuse:
    std: standard deviation
"""

def calc_std_3(values):
    n = len(values)
    values = values - np.mean(values) #subtract mean to make the mean of the values 0
    s = np.sum(values**2)
    
    inside_sqrt = (1/(n-1)) * (s - n * (np.mean(values))**2)
    
    if inside_sqrt < 0:
        print ('Error: negative value inside the square root')
    
    std = np.sqrt(abs(inside_sqrt))
    
    return std

print('\n', 'Question 1.d')
rel_error_3_1e7 = relative_error(vals2, calc_std_3)
print('The relative error using the modified eq. (2) for mean = 1e7 is ', rel_error_3_1e7)
    
    
    
