# -*- coding: utf-8 -*-
"""
Created on Sun Mar  2 10:30:46 2025

@author: isabe
"""

import numpy as np
import matplotlib.pyplot as plt
from ddeint import ddeint
from scipy.optimize import fsolve
from scipy.optimize import root

#parameters
c1, c2, c3, c4 = 10, -10, 10, 2.0
P, Q = -2.5, -4.0
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

plt.figure(figsize=(14,9))
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

def steady_state(vars):
    u1, v1, u2, v2 = vars

    eq1 = u1 - f(c1 * u1 + c2 * v1 + P + eps * u2)
    eq2 = v1 - f(c3 * u1 + c4 * v1 + Q)
    eq3 = u2 - f(c1 * u2 + c2 * v2 + P + eps * u1)
    eq4 = v2 - f(c3 * u2 + c4 * v2 + Q)

    return [eq1, eq2, eq3, eq4]

initial_guess = [0.5, 0.5, 0.5, 0.5]

# Solve for steady states
steady_states = fsolve(steady_state, initial_guess)

u1_star, v1_star, u2_star, v2_star = steady_states

# Display results
print(f"Steady state values:")
print(f"u1* = {u1_star:.4f}")
print(f"v1* = {v1_star:.4f}")
print(f"u2* = {u2_star:.4f}")
print(f"v2* = {v2_star:.4f}")

def build_Jacobian(λ, u1_star, v1_star, u2_star, v2_star):
    phi1 = f_prime(c1 * u1_star + c2 * v1_star + P + eps * u2_star)
    phi2 = f_prime(c3 * u1_star + c4 * v1_star + Q)
    phi3 = f_prime(c1 * u2_star + c2 * v2_star + P + eps * u1_star)
    phi4 = f_prime(c3 * u2_star + c4 * v2_star + Q)

    J = np.array([
        [-1 + c1 * phi1 * np.exp(-λ * τ1), c2 * phi1 * np.exp(-λ * τ2), eps * phi1 * np.exp(-λ * τ3), 0],
        [c3 * phi2 * np.exp(-λ * τ2), -1 + c4 * phi2 * np.exp(-λ * τ1), 0, 0],
        [eps * phi3 * np.exp(-λ * τ3), 0, -1 + c1 * phi3 * np.exp(-λ * τ1), c2 * phi3 * np.exp(-λ * τ2)],
        [0, 0, c3 * phi4 * np.exp(-λ * τ2), -1 + c4 * phi4 * np.exp(-λ * τ1)]
    ])
    
    return J

def characteristic_equation(lmbda_complex, u1_star, v1_star, u2_star, v2_star):
    λ = lmbda_complex[0] + 1j * lmbda_complex[1]
    J = build_Jacobian(λ, u1_star, v1_star, u2_star, v2_star)
    I = np.eye(4, dtype=complex)
    char_eq = np.linalg.det(J - λ * I)
    return [char_eq.real, char_eq.imag]

# Initial guess for eigenvalue
lambda_guess = [0, 0]

# Solve for λ
result = root(characteristic_equation, lambda_guess, args=(u1_star, v1_star, u2_star, v2_star))

# Extract λ = real + imag
λ = result.x[0] + 1j * result.x[1]
print(f"Computed eigenvalue: λ = {λ.real:.4f} + {λ.imag:.4f}j")
