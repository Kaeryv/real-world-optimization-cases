import numpy as np
import matplotlib.pyplot as plt

from .misc import array_factorcir, trapezoidalcir
def display_plot(gbest, phi_desired, distance):
    dim = len(gbest)
    pi = np.pi
    phi = np.linspace(0, 360, 5000)
    yax = np.zeros(5000)

    yax[0] = array_factorcir(gbest, (pi/180)*phi[0], phi_desired, distance, dim)
    maxi = yax[0]

    for i in range(1, 5000):
        yax[i] = array_factorcir(gbest, (pi/180)*phi[i], phi_desired, distance, dim)
        if maxi < yax[i]:
            maxi = yax[i]

    for i in range(5000):
        yax[i] = yax[i] / maxi
        yax[i] = 20 * np.log10(yax[i])

    plt.plot(phi, yax, 'g')
    plt.xlabel('Azimuth angle (deg)')
    plt.ylabel('Gain (db)')
    plt.show()

def antennafunccircular(x1, null, phi_desired, distance):
    global directivity

    pi = 3.141592654
    dim = len(x1)
    y = 0

    num_null = len(null)
    num1 = 300
    phi = np.linspace(0, 360, num1)
    phizero = 0
    yax = np.zeros(num1)

    yax[0] = array_factorcir(x1, (pi/180)*phi[0], phi_desired, distance, dim)
    maxi = yax[0]
    phi_ref = 1

    for i in range(1, num1):
        yax[i] = array_factorcir(x1, (pi/180)*phi[i], phi_desired, distance, dim)
        if maxi < yax[i]:
            maxi = yax[i]
            phizero = phi[i]
            phi_ref = i + 1

    sidelobes = np.zeros(num1)
    sllphi = np.zeros(num1)
    count = 0

    if yax[0] > yax[num1-1] and yax[0] > yax[1]:
        count += 1
        sidelobes[count] = yax[0]
        sllphi[count] = phi[0]

    if yax[num1-1] > yax[0] and yax[num1-1] > yax[num1-2]:
        count += 1
        sidelobes[count] = yax[num1-1]
        sllphi[count] = phi[num1-1]

    for i in range(1, num1-1):
        if yax[i] > yax[i+1] and yax[i] > yax[i-1]:
            count += 1
            sidelobes[count] = yax[i]
            sllphi[count] = phi[i]

    sidelobes = np.sort(sidelobes)[::-1]

    upper_bound = 180
    lower_bound = 180

    y = sidelobes[1] / maxi
    sllreturn = 20 * np.log10(y)

    for i in range(1, num1//2 + 1):
        if phi_ref + i > num1 - 1:
            upper_bound = 180
            break
        tem = yax[phi_ref + i - 1]
        if yax[phi_ref + i - 1] < yax[phi_ref + i - 2] and yax[phi_ref + i - 1] < yax[phi_ref + i]:
            upper_bound = phi[phi_ref + i - 1] - phi[phi_ref - 1]
            break

    for i in range(1, num1//2 + 1):
        if phi_ref - i < 2:
            lower_bound = 180
            break
        tem = yax[phi_ref - i - 1]
        if yax[phi_ref - i - 1] < yax[phi_ref - i - 2] and yax[phi_ref - i - 1] < yax[phi_ref - i]:
            lower_bound = phi[phi_ref - 1] - phi[phi_ref - i - 1]
            break

    bwfn = upper_bound + lower_bound

    y1 = 0
    for i in range(num_null):
        y1 += array_factorcir(x1, null[i], phi_desired, distance, dim) / maxi

    y2 = 0
    uavg = trapezoidalcir(x1, 0, 2 * pi, 50, phi_desired, distance, dim)
    y2 = np.abs(2 * pi * maxi**2 / uavg)

    directivity = 10 * np.log10(y2)

    y3 = np.abs(phizero - phi_desired)
    if y3 < 5:
        y3 = 0

    y = 0
    if bwfn > 80:
        y += np.abs(bwfn - 80)

    y = sllreturn + y + y1 + y3

    return y, sllreturn, bwfn
