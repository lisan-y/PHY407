"""
PHY407 Lab 02

Functions for Lab 2; used in Lab02_Q2a.py, Lab02_Q2b.py, Lab02_Q3.py, and Q2 c (code not submitted)

Author: Lisa Nasu-Yu, 24 Sept 2021
"""


def trapezoid(f, a, b, N):
    """
    Evalutes integrals using trapezoid rule with N slices using the formula
    I(a,b) = h * [1/2 * f(a) + 1/2 * f(b) + sum(f(a + kh)] for k in [1, N-1],
    derived from manipulation of the Taylor expansion of f.
    :param f: [NumPy function] integrand
    :param a: [float] lower bound of integration
    :param b: [float] upper bound of integration
    :param N: [int] number of slices
    :return: [float] solution to integral with trapezoid rule for N slices
    """
    h = (b - a) / N  # width of slices
    s = 0.5 * f(a) + 0.5 * f(b)  # end slices
    for k in range(1, N):  # add each interior slice
        s += f(a + k * h)
    return h * s


def simpsons(f, a, b, N):
    """
    Evalutes integrals using Simpson's rule with N slices using the formula
    I(a,b) = 1/3 * h * [f(a) + f(b) + 4 * sum(f(a + k*h)) + 2 * sum(f(a + m*h))]
    for odd k in [1, N-1] and even m in [2, N-2],
    derived from manipulation of the Taylor expansion of f.
    :param f: [NumPy function] integrand
    :param a: [float] lower bound of integration
    :param b: [float] upper bound of integration
    :param N: [int] number of slices
    :return: [float] solution to integral with Simpson's rule for N slices
        """
    h = (b - a) / N  # width of slices
    s = f(a) + f(b)  # end slices
    for k in range(1, N, 2):  # odd terms
        s += 4 * f(a + k * h)
    for k in range(2, N, 2):  # even terms
        s += 2 * f(a + k * h)
    return (h / 3) * s
