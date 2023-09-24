import numpy as np
from numpy.linalg import solve

A = np.array([[3, 5, 6, 2], [2, 0, 1, 3], [4, 1, 5, 6], [-5, 3, 5, 3]], float)
v = np.array([24, 7, 13, -7], float)

x = solve(A, v)
print(x)

print(np.linalg.eigvalsh(A))

n=1
while 30/(2**n) > 1e-4:
    n += 1

print(n, 30/(2**(n-1)), 30/(2**n))
