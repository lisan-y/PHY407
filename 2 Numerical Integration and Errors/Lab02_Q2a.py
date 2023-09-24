"""
PHY407 Lab 02
Q2 a
ii) Computes the integral of 4/(1+x**2) using the trapezoid and Simpson's rules.
iii) Determines the number of slices necessary to give an error of order -9 for each integration method,
and measures the computation time of the corresponding integration.
iv) Estimates the error for the trapezoidal rule with 32 slices using the practical estimation of errors.

Author: Lisa Nasu-Yu, 24 Sept 2021
"""

from time import time
import numpy as np
from functions_lab02 import *  # make sure in same folder


# Define integrand, and its derivative
def f(x):
    """
    Evaluate the function 4 / (1 + x**2) at x
    :param x: [float] value at which to evaluate function
    :return: [float] value of function at x
    """
    return 4 / (1 + x ** 2)


def df(x):
    """
    Evaluate the derivative of the function 4 / (1 + x**2) at x
    :param x: [float] value at which to evaluate the derivative
    :return: [float] value of the derivative x
    """
    return -8 * x / (1 + x ** 2) ** 2


print("Q2a ii")

# Define boundaries and slices
N = 4  # number of slices
a = 0.0  # Lower bound of integral
b = 1.0  # upper bound of integral

print("The value with trapezoidal rule and {} slices is {}".format(N, trapezoid(f, a, b, N)))
print("The value with Simpson's rule and {} slices is {}".format(N, simpsons(f, a, b, N)))


print("\nQ2a iii")

# Trapezoidal Rule
# Find N at which trapezoidal rule has error of order -9
# Define number of slices as 2 ** (n = 2)
n = 2
h = (b - a) / (2 ** n)  # slice width
# Calculate error with Euler-Maclaurin formula
error = (1 / 12) * h ** 2 * (df(a) - df(b))
# calculate error for increasing n until the error is greater than 1e-8
while error > 1e-8:
    n += 1
    h = (b - a) / (2 ** n)  # slice width
    error = (1 / 12) * h ** 2 * (df(a) - df(b))

# time 10 integrations with error O(-9)
start = time()
for i in range(100):
    trapezoid(f, a, b, 2 ** n)
end = time()
diff = end - start
# divide by 100 to get average time of a single integration
diff /= 100

print("2 ** {} slices gives an error of order -9 for the trapezoidal rule."
      "\n\tThe error for 2 ** {} slices is {:.3g}".format(n, n, error))
print("\tIt takes an average of {:.3g}s to compute the integral for 2 ** {} slices".format(diff, n))

# Simpson's Rule
# Find N at which Simpson's rule has error of order -9, using same method as above
n = 2
# Calculate error as the difference between Simpson's rule value from actual value, pi
error = np.abs(np.pi - simpsons(f, a, b, 2 ** n))
while error > 1e-8:
    n += 1
    error = np.abs(np.pi - simpsons(f, a, b, 2 ** n))

# time 10 integrations with error O(-9)
start = time()
for i in range(100):
    simpsons(f, a, b, 2 ** n)
end = time()
diff = end - start
# divide by 100 to get average time of a single integration
diff /= 100

print("2 ** {} slices gives an error of order -9 for the Simpson's rule."
      "\n\tThe error for 2 ** {} slices is {:.3g}".format(n, n, error))
print("\tIt takes an average of {:.3g}s to compute the integral for 2 ** {} slices".format(diff, n))

print("\nQ2a iv")
# Define N1 and N2, number of steps
N1 = 16
N2 = 32

# Evaluate integral with trapezoid rule for N1 and N2 steps
I1 = trapezoid(f, a, b, N1)
I2 = trapezoid(f, a, b, N2)

# Evaluate error for N2 slices with practical estimation of errors method
e2 = (I2 - I1) / 3

print("The error estimation for the trapezoidal rule with N2 = 32 slices is {:.3g}".format(e2))
