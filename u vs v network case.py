# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 15:24:08 2025

@author: isabe
"""

import numpy as np
import matplotlib.pyplot as plt
from ddeint import ddeint
from scipy.optimize import fsolve
from scipy.optimize import root

#parameters
c1, c2, c3, c4 = 10, -10, 10, 2.0
P, Q = -2.6, -4.0
τ1, τ2, τ3 = 0.5, 0.5, 2.0
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
    u1_tau1, u2_tau1 = Y(t - τ1)[0], Y(t - τ1)[2]
    v1_tau1, v2_tau1 = Y(t - τ1)[1], Y(t - τ1)[3]
    u1_tau2, u2_tau2 = Y(t - τ2)[0], Y(t - τ2)[2]
    v1_tau2, v2_tau2 = Y(t - τ2)[1], Y(t - τ2)[3]
    u2_tau3 = Y(t - τ3)[2]
    u1_tau3 = Y(t - τ3)[0]

    #mass 1
    du1dt = -u1 + f(c1 * u1_tau1 + c2 * v1_tau2 + P + eps * u2_tau3)
    dv1dt = -v1 + f(c3 * u1_tau2 + c4 * v1_tau1 + Q)

    #mass 2
    du2dt = -u2 + f(c1 * u2_tau1 + c2 * v2_tau2 + P + eps * u1_tau3)
    dv2dt = -v2 + f(c3 * u2_tau2 + c4 * v2_tau1 + Q)

    return np.array([du1dt, dv1dt, du2dt, dv2dt])

t = np.linspace(0, 50, 1000)
result = ddeint(system, hist_func, t)

#extract solutions
u1, v1, u2, v2 = result.T

plt.figure(figsize=(14, 10))
plt.subplot(2, 1, 1)
plt.plot(u1, v1, label="u1 vs v1 (Mass 1)")
plt.legend(fontsize=20)
plt.grid()

plt.subplot(2, 1, 2)
plt.plot(u2, v2, label="u2 vs v2 (Mass 2)")
plt.legend(fontsize=20)
plt.grid()

plt.xlabel("u(t)", fontsize=18)
plt.ylabel("v(t)", fontsize=18)

plt.show()

