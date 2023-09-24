#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 10:30:35 2021

@author: emmajarvis
"""
import numpy as np
from scipy.constants import h, c, k

def f(x):
    return 5*np.exp(-x)+x-5
def df(x):
    return -5*np.exp(-x)+1

def binary_search(f, x1, x2):
    err = x2-x1
    i = 0

    while err >= 1e-6:
        f1 = f(x1)
        f2 = f(x2)
        x_mid = 1/2*(x2+x1)
        if f1*f(x_mid)>0:
            x1 = x_mid
        elif f2*f(x_mid)>0:
            x2 = x_mid
        err = x2-x1
        i = i + 1
    print(i)
     
    return x_mid, err

def relaxation(f, x1):
    x = x1
    i = 0
    err = f(x)
    while abs(err) >= 1e-6:
        err = f(x)
        x = 5-5*np.exp(-x)
        i = i + 1
    print(i)
    return x, err


def newtons(f, x1):
    x = x1
    err = f(x)
    i = 0
    while abs(err) >= 1e-6:
        err = f(x)
        x = x - f(x)/df(x)
        i = i + 1
    print(i)
    return x, err



x, err = binary_search(f, 4, 6)
print('x: ', x, 'error: ', err)

def displacement_constant(x):
    return h*c/(k*x)

b = displacement_constant(x)

print('The displacement constant is', b, 'Km')

lam = 500e-9 #m
T = b/lam

print('The temperature of the sun is ', T, 'K')

start = 6

print('Starting Index')


print('\n Binary Search')
x_b, err = binary_search(f, 4, start)
print('x: ', x_b, 'error:', err)

print('\n Relaxation Method')
x_r, err = relaxation(f, start)
print('x: ', x_r, 'error:', err)

print('\n Newtons Method')
x_n, err = newtons(f, start)
print('x: ', x_n, 'error:', err)