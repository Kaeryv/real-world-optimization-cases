import numpy as np
from math import sin, cos, exp, pi, floor
from scipy.io import loadmat

def function_selector(x, num=1):
    return globals()[f"rw{num}"](x)

def get_boundaries(num):
    bnds = None
    if num == 1:
        bnds = np.zeros((2, 6))
        bnds[0, :] = -6.4
        bnds[1, :] = +6.35
    elif num == 2 or num == 5 or num == 6:
        bnds = list()
        bnds.append((0, 4))
        bnds.append((0, 4))
        bnds.append((0, pi))

        for i in range(4, 31):
            alpha = floor((i-4)/3)
            bnds.append((-4 - alpha/4, 4+alpha/4))
        bnds = np.asarray(bnds).T
    elif num == 3:
        bnds = np.array([[0.6], [0.9]])
    elif num == 4:
        bnds = np.array([[0], [5]])

    elif num == 7:
        bnds = np.zeros((2, 20))
        bnds[0, :] = 0
        bnds[1, :] = 2*pi
        

    else:
        print("Unsupported problem.")

    return bnds

def diffsolv(x, t, u):
    c = loadmat('c_bifunc_data.mat')["c"]  # Assuming the data file is in the same directory
    ml = np.array([1, u[0], u[0]**2, u[0]**3])
    mlt = np.tile(ml, (10, 1))
    k = np.sum(c * mlt, axis=1)

    dy = np.zeros(7)
    dy[0] = -k[0] * x[0]
    dy[1] = k[0] * x[0] - (k[1] + k[2]) * x[1] + k[3] * x[4]
    dy[2] = k[1] * x[1]
    dy[3] = -k[5] * x[3] + k[4] * x[4]
    dy[4] = k[2] * x[1] + k[5] * x[3] - (k[3] + k[4] + k[7] + k[8]) * x[4] + k[6] * x[5] + k[9] * x[6]
    dy[5] = k[7] * x[4] - k[6] * x[5]
    dy[6] = k[8] * x[4] - k[9] * x[6]

    return dy


def intgrl(x, t, u):
    dy = np.zeros(2)
    dy[0] = -(2 + u) * (x[0] + 0.25) + (x[1] + 0.5) * exp(25 * x[0] / (x[0] + 2))
    dy[1] = 0.5 - x[1] - (x[1] + 0.5) * exp(25 * x[0] / (x[0] + 2))
    return dy

def rw1(x):
    '''
        FM Sound wave
    '''
    d = len(x)  # Assuming x is defined elsewhere

    if len(x) != 6:
        print('dimension-size should be six.')
    
    theta = 2 * np.pi / 100
    f = 0
    
    for t in range(101):
        y_t = x[0] * sin(x[1] * t * theta + x[2] * sin(x[3] * t * theta + x[4] * sin(x[5] * t * theta)))
        y_0_t = 1 * sin(5 * t * theta - 1.5 * sin(4.8 * t * theta + 2 * sin(4.9 * t * theta)))
        f += (y_t - y_0_t)**2

    return f


def rw2(x):
    '''
        Lennard Jones potential
    '''
    x = np.asarray(x) if isinstance(x, list) else x
    if x.shape[0] % 3 != 0:
        print('x passed to this function must be n dimensional array where, n is perfectly divisible by 3.')
        return np.nan
    n = x.shape[0] // 3
    r = np.zeros((n, n))
    x = np.reshape(x, (n, 3))
    v = 0
    a = np.ones((n, n))
    b = 2*np.ones((n, n))
    for i in range(n-1):
        for j in range(i+1, n):
            r[i,j] = np.linalg.norm(x[i, :] - x[j, :])
            v += (a[i, j] / r[i,j]**12 - b[i, j] / r[i,j]**6)

    return v

from scipy.integrate import odeint

def rw3(x):
    '''
        BiF Catalyst blend.
    '''
    x = np.asarray(x) if isinstance(x, list) else x
    tol = 1.0e-01
    tspan = [0, 0.78]
    tspan = np.linspace(tspan[0], tspan[1], 100)
    yo = [1, 0, 0, 0, 0, 0, 0]

    # Solving ODEs
    sol = odeint(diffsolv, yo, tspan, args=(x,), rtol=tol)

    # Extracting the final value
    f = sol[-1, -1] * 1e3

    return f


def rw4(x):
    '''
        Stirred tank reactor
    '''
    x = np.asarray(x) if isinstance(x, list) else x
    tol = 1.0e-01
    tspan = [0, 0.78]
    yo = np.array([0.09, 0.09])
    u = x  # u should be passed here.
    options = {'atol': tol, 'rtol': tol}
    
    Y = odeint(intgrl, yo, np.linspace(tspan[0], tspan[1], 100), hmax=0.01, args=(x,), **options)
    f = np.sum(np.sum(Y**2, axis=1) + 0.1 * u**2)

    return f

def rw5(x):
    '''
        Tersoff Si (B)
    '''
    x = np.asarray(x) if isinstance(x, list) else x
    p = x.shape
    if p[0] % 3 != 0:
        print('x passed to this function must be n dimensional array where, n is perfectly divisible by 3.')

    NP = p[0] // 3
    x = np.reshape(x, (NP, 3))
    R1 = 3.0
    R2 = 0.2
    A = 3.2647e+3
    B = 9.5373e+1
    lemda1 = 3.2394
    lemda2 = 1.3258
    lemda3 = 1.3258
    c = 4.8381
    d = 2.0417
    n1 = 22.956
    gama = 0.33675
    h = 0
    E = np.zeros(NP)
    r = np.zeros((NP, NP))
    fcr = np.zeros((NP, NP))
    VRr = np.zeros((NP, NP))
    VAr = np.zeros((NP, NP))

    for i in range(NP):
        for j in range(NP):
            r[i, j] = np.sqrt(np.sum((x[i, :] - x[j, :])**2))
            if r[i, j] < (R1 - R2):
                fcr[i, j] = 1
            elif r[i, j] > (R1 + R2):
                fcr[i, j] = 0
            else:
                fcr[i, j] = 0.5 - 0.5 * sin(np.pi / 2 * (r[i, j] - R1) / R2)

            VRr[i, j] = A * exp(-lemda1 * r[i, j])
            VAr[i, j] = B * exp(-lemda2 * r[i, j])

    for i in range(NP):
        for j in range(NP):
            if i == j:
                continue

            jeta = np.zeros((NP, NP))
            for k in range(NP):
                if i == k or j == k:
                    continue

                rd1 = np.linalg.norm(x[i, :] - x[k, :], axis=0)
                rd3 = np.linalg.norm(x[k, :] - x[j, :], axis=0)
                rd2 = np.linalg.norm(x[i, :] - x[j, :], axis=0)
                ctheta_ijk = (rd1**2 + rd2**2 - rd3**3) / (2 * rd1 * rd2)
                G_th_ijk = 1 + (c**2) / (d**2) - (c**2) / (d**2 + (h - ctheta_ijk)**2)
                jeta[i, j] += fcr[i, k] * G_th_ijk * np.exp(lemda3**3 * (r[i, j] - r[i, k])**3)

            Bij = (1 + (gama * jeta[i, j])**n1)**(-0.5 / n1)
            E[i] += fcr[i, j] * (VRr[i, j] - Bij * VAr[i, j]) / 2

    return np.nansum(E)



def rw6(x):
    '''
        Tersoff potential for Si (C)
    '''
    p = x.shape
    if p[0] % 3 != 0:
        print('x passed to this function must be n dimensional array where, n is perfectly divisible by 3.')

    NP = p[0] // 3
    x = np.reshape(x, (NP, 3))
    R1 = 2.85
    R2 = 0.15
    A = 1.8308e+3
    B = 4.7118e+2
    lemda1 = 2.4799
    lemda2 = 1.7322
    lemda3 = 1.7322
    c = 1.0039e+05
    d = 1.6218e+01
    n1 = 7.8734e-01
    gama = 1.0999e-06
    h = -5.9826e-01
    E = np.zeros(NP)
    r = np.zeros((NP, NP))
    fcr = np.zeros((NP, NP))
    VRr = np.zeros((NP, NP))
    VAr = np.zeros((NP, NP))

    for i in range(NP):
        for j in range(NP):
            r[i, j] = np.sqrt(np.sum((x[i, :] - x[j, :])**2))
            if r[i, j] < (R1 - R2):
                fcr[i, j] = 1
            elif r[i, j] > (R1 + R2):
                fcr[i, j] = 0
            else:
                fcr[i, j] = 0.5 - 0.5 * np.sin(np.pi / 2 * (r[i, j] - R1) / R2)

            VRr[i, j] = A * np.exp(-lemda1 * r[i, j])
            VAr[i, j] = B * np.exp(-lemda2 * r[i, j])

    for i in range(NP):
        for j in range(NP):
            if i == j:
                continue

            jeta = np.zeros((NP, NP))
            for k in range(NP):
                if i == k or j == k:
                    continue

                rd1 = np.sqrt(np.sum((x[i, :] - x[k, :])**2))
                rd3 = np.sqrt(np.sum((x[k, :] - x[j, :])**2))
                rd2 = np.sqrt(np.sum((x[i, :] - x[j, :])**2))
                ctheta_ijk = (rd1**2 + rd2**2 - rd3**3) / (2 * rd1 * rd2)
                G_th_ijk = 1 + (c**2) / (d**2) - (c**2) / (d**2 + (h - ctheta_ijk)**2)
                jeta[i, j] += fcr[i, k] * G_th_ijk * np.exp(lemda3**3 * (r[i, j] - r[i, k])**3)

            Bij = (1 + (gama * jeta[i, j])**n1)**(-0.5 / n1)
            E[i] += fcr[i, j] * (VRr[i, j] - Bij * VAr[i, j]) / 2

    return np.nansum(E)


def rw7(x):
    '''
        Sprd Spectrum rad Phase
    '''
    d = len(x)
    hsum = np.zeros((2*d-1)*2)
    var = 2 * d - 1

    for kk in range(1, 2 * var + 1):
        if kk % 2 != 0:
            i = (kk + 1) // 2
            for j in range(i, d + 1):  # fi(2i-1)X
                summ = 0
                for i1 in range(abs(2 * i - j - 1) + 1, j + 1):
                    summ += x[i1 - 1]
                hsum[kk-1] += np.cos(summ)
        else:
            i = kk // 2
            for j in range(i + 1, d + 1):  # fi(2i)X
                summ = 0
                for i1 in range(abs(2 * i - j) + 1, j + 1):
                    summ = x[i1 - 1] + summ
                hsum[kk-1] += np.cos(summ)
            hsum[kk-1] += 0.5

    return np.max(hsum)



# def rw8(x):
#     # Define the function for Transmission Network Expansion Planning
#     sw = np.ceil(x)
#     # Assuming Candidate and Linedata are defined elsewhere
#     n1 = len(Linedata[:, 1])
#     sw1 = sw.astype(int)

#     for k in range(len(sw1)):
#         Linedata[n1 + k, :] = Candidate[sw1[k], :]

#     n_orginalLine = n1
#     n = len(Pgen)
#     B = np.zeros((n, n))
#     Nline = len(Linedata[:, 1])
#     Xline = Linedata[:, 4]
#     pijmax = Linedata[:, 6]
#     Tap = np.ones((n, n))

#     for C in range(Nline):
#         bline = 1 / Xline[C]
#         k = Linedata[C, 2]
#         m = Linedata[C, 3]
#         B[k, m] -= bline
#         B[m, k] = B[k, m]
#         B[k, k] += bline
#         B[m, m] += bline

#     B[0, 0] = 10000000
#     X = np.linalg.inv(B)
#     delP = Pgen - Pload
#     delP = delP.T
#     delta = np.dot(X, delP)
#     pij = np.zeros(Nline)

#     for k in range(Nline):
#         i = Linedata[k, 2]
#         j = Linedata[k, 3]
#         pij[k] = (delta[i] - delta[j]) / Xline[k]

#     PIPbase = 0.0
#     f = np.sum(Linedata[n_orginalLine + 1:, 7]) + 30
#     pen = 0

#     for i in range(len(Linedata[:, 1])):
#         pen += 5000 * max((np.abs(pij[i]) - Linedata[i, 6]), 0)

#     for i in range(len(Candidate[:, 1])):
#         a = np.where(sw == i)[0]
#         if len(a) > 3:
#             pen += 1000

#     f += pen

#     return f




