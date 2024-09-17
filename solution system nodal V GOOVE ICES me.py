# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 14:07:28 2024

@author: Admin
"""

import numpy as np

# Constants
sigma = 5.67e-8  # Stefan-Boltzmann constant
T1 = 245  # Temperature in Kelvin
Ts = 4  # Temperature in Kelvin (assuming this is for space)

# Emissivity and surface areas
epsilon = 0.023
epsilon_out = 0.8

A1 = 4  # Surface A1
A2 = 4  # Surface A2
A3 = 4  # Surface A3

# View factors
F1s = 0.052  # View factor from surface 1 to space
F2s = 0.052  # View factor from surface 2 to space
F3s = 0.8    # View factor from surface 3 to space
F12 = 0.948  # View factor from surface 1 to surface 2
F23 = 0.948  # View factor from surface 2 to surface 3

# Thermal resistances
R1s = 1 / (F1s * A1)
R2s = 1 / (F2s * A2)
R3s = 1 / (F3s * A3)

R12 = 1 / (F12 * A1)
R23 = 1 / (F23 * A2) 

R1 = (1 - epsilon) / (A1 * epsilon)
R2 = (1 - epsilon) / (A2 * epsilon)
R3 = (1 - epsilon_out) / (A3 * epsilon_out)

# Matrix A
A = np.array([
    [(-1/R1 - 1/R12 - 1/R1s), 1/R12, 0, 1/R1s, 1/R1, 0, 0, 0],
    [1/R12, (-1/R2 - 1/R23 - 1/R12 - 1/R2s), 1/R23, 1/R2s, 0, 1/R2, 0, 0],
    [0, 1/R23, (-1/R3 - 1/R23 - 1/R3s), 1/R3s, 0, 0, 1/R3, 0],
    [-1/R1s, -1/R2s, -1/R3s, (1/R1s + 1/R2s + 1/R3s), 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 1, 0, 0, 0, -1],
    # [1, 0, 0, 0, -1, 0, 0, 0]
    [1/R1, 1/R2, 1/R3, 0, -1/R1, -1/R2, -1/R3, 0]
])

# Vector Y
Y = np.array([
    0,
    0,
    0,
    0,
    sigma * T1**4,
    sigma * Ts**4,
    0,
    0
])

# Solve the system of equations
X = np.linalg.solve(A, Y)

print(X)

T1 = (X[0] / sigma)**(1/4)
T2 = (X[1] / sigma)**(1/4)
T3 = (X[2] / sigma)**(1/4)
T4 = (X[3] / sigma)**(1/4)
T5 = (X[4] / sigma)**(1/4)
T6 = (X[5] / sigma)**(1/4)
T7 = (X[6] / sigma)**(1/4)
T8 = (X[7] / sigma)**(1/4)


print(f"Temperature T1: {T1:.2f} K")
print(f"Temperature T2: {T2:.2f} K")
print(f"Temperature T3: {T3:.2f} K")
print(f"Temperature Ts: {T4:.2f} K")
print(f"Temperature T1: {T5:.2f} K")
print(f"Temperature T2: {T6:.2f} K")
print(f"Temperature T3: {T7:.2f} K")
print(f"Temperature Ts: {T8:.2f} K")
