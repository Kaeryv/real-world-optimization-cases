import numpy as np
from numba import njit

@njit
def array_factorcir(x1, phi, phi_desired, distance, dim):
    pi = 3.141592654
    y = 0
    y1 = 0

    for i1 in range(1, dim//2 + 1):
        delphi = 2 * pi * (i1 - 1) / dim
        shi = np.cos(phi - delphi) - np.cos(phi_desired * (pi / 180) - delphi)
        shi = shi * dim * distance
        y += x1[i1 - 1] * np.cos(shi + x1[dim//2 + i1 - 1] * (pi / 180))

    for i1 in range(dim//2 + 1, dim + 1):
        delphi = 2 * pi * (i1 - 1) / dim
        shi = np.cos(phi - delphi) - np.cos(phi_desired * (pi / 180) - delphi)
        shi = shi * dim * distance
        y += x1[i1 - dim//2 - 1] * np.cos(shi - x1[i1 - 1] * (pi / 180))

    for i1 in range(1, dim//2 + 1):
        delphi = 2 * pi * (i1 - 1) / dim
        shi = np.cos(phi - delphi) - np.cos(phi_desired * (pi / 180) - delphi)
        shi = shi * dim * distance
        y1 += x1[i1 - 1] * np.sin(shi + x1[dim//2 + i1 - 1] * (pi / 180))

    for i1 in range(dim//2 + 1, dim + 1):
        delphi = 2 * pi * (i1 - 1) / dim
        shi = np.cos(phi - delphi) - np.cos(phi_desired * (pi / 180) - delphi)
        shi = shi * dim * distance
        y1 += x1[i1 - dim//2 - 1] * np.sin(shi - x1[i1 - 1] * (pi / 180))

    y = y * y + y1 * y1
    y = np.sqrt(y)

    return y


import numpy as np
def trapezoidalcir(x2, upper, lower, N1, phi_desired, distance, dim):
    # This function performs integration by trapezoidal rule
    h = (upper - lower) / N1
    x1 = lower
    y = np.abs(np.real(array_factorcir(x2, lower, phi_desired, distance, dim))**2 * np.sin(lower - np.pi/2))

    for i in range(2, N1 + 2):
        x1 += h
        y = np.append(y, np.abs(np.real(array_factorcir(x2, x1, phi_desired, distance, dim))**2 * np.sin(x1 - np.pi/2)))

    s = 0
    for i in range(1, N1 + 2):
        if i == 1 or i == N1 + 2:
            s += y[i - 1]
        else:
            s += 2 * y[i - 1]

    q = (h / 2) * s
    return q
