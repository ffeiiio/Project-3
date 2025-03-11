# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 09:42:22 2025

@author: isabe
"""

import numpy as np
import matplotlib.pyplot as plt
from ddeint import ddeint
from scipy.optimize import fsolve
from scipy.optimize import root

#parameters
c1, c2, c3, c4 = 10, -10, 10, 2.0
P, Q = -2.4, -4.0
eps = 1.0
ß = 1.0


def f(z):
    return 1 / (1 + np.exp(-ß * z))

def f_prime(z, ß=1):
    fz = f(z)
    return ß * fz * (1 - fz)

#history function
def hist_func(t):
    return np.array([0.1, 0.1, 0, 0])

#system of DDEs
def system(Y, t):
    u1, v1, u2, v2 = Y(t)

    #mass 1
    du1dt = -u1 + f(c1 * u1 + c2 * v1 + P + eps * u2)
    dv1dt = -v1 + f(c3 * u1 + c4 * v1 + Q)

    #mass 2
    du2dt = -u2 + f(c1 * u2 + c2 * v2 + P + eps * u1)
    dv2dt = -v2 + f(c3 * u2 + c4 * v2 + Q)

    return np.array([du1dt, dv1dt, du2dt, dv2dt])

t = np.linspace(0, 50, 1000)
result = ddeint(system, hist_func, t)

#extract solutions
u1, v1, u2, v2 = result.T

plt.figure(figsize=(14, 9))
plt.subplot(2, 1, 1)
plt.plot(t, u1, label="u1")
plt.plot(t, v1, label="v1")
plt.legend(fontsize=20)
plt.grid()

plt.subplot(2, 1, 2)
plt.plot(t, u2, label="u2")
plt.plot(t, v2, label="v2")
plt.legend(fontsize=20)
plt.grid()

plt.xlabel("Time", fontsize=18)

plt.show()