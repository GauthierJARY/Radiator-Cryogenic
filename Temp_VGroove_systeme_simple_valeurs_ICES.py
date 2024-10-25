# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 10:36:59 2024

@author: jaryg

Code réalisé dans le cadre d'un stage au CNES. 
Il s'agit d'estimer les facteurs de Gebhart dans une configuration de radiateurs cryogéniques V-Groove (cas applicatif du satellite Planck par exemple)
Pour cela, la matrice et la résolution analytique des facteurs de vue a été réalisée, et ce afin de comparer avec les résultats donnés par une autre méthode, comme le tir de rayons.
"""

import numpy as np
import math as m

# Constants and parameters

sigma = 5.8e-8 # W/m2, Stefan-Boltzmann constant
eps1 = 0.023 # emissivity of V Groove panel 1
eps2 = 0.023 # emissivity of V Groove panel 2
eps3 = 0.023 # emissivity of V Groove panel 3 (inner side)
eps3f = 0.8 # emissivity of V Groove panel 3 (outer side)
eps4 = 1 # emissivity of the cryostat (approximated as space)
theta = 6 # degrees
theta = theta/360 * 2 * m.pi # convert degrees to radians
theta1 = theta
theta2 = theta
theta3 = theta
A1 = 0.0268973 # m2, area of V Groove panel 1
A2 = 0.0255773 * 2 # m2, area of V Groove panel 2 (doubled)
A3 = 0.0244109 * 2 # m2, area of V Groove panel 3 (doubled)

T1 = 245 # K, temperature of V Groove panel 1
Ts = 0 # K, temperature of the cryostat/space

# View factor / Form factor
# Note: this is a purely geometric reasoning, without considering ray rebounds

# view from panel 1
F12 = 1 - m.sin(theta1/2)
F14 = m.sin(theta1/2)
# view from panel 2
F21 = A1/A2 * F12
F24 = 2 * m.sin(theta2/2)
F23 = 1 - F24 - F21
# view from panel 3
F32 = F23 * A2/A3 # Note: calculated using the reflectivity formula; adjusting values can improve power ratios
F34 = 1 - F32 
# view from the cryostat / space
F41 = 0.0 # negligible, according to the reciprocity formula due to space being nearly infinite
F42 = 0.0
F43 = 0.0
F3o4 = 1
# internal and external view of panel 3
F3i2 = 1 - m.sin(theta/2)
F23i = 1 - m.sin(theta/2)

# Note: If the view factors are incorrect or not normalized, they need to be adjusted. Normalization follows:
norm1 = F12 + F14
F12 = F12 / norm1
F14 = F14 / norm1
# norm2 = F21 + F23 + F24
# F21 = F21 / norm2
# F23 = F23 / norm2
# F24 = F24 / norm2
norm3 = F32 + F34
F32 = F32 / norm3
F34 = F34 / norm3
norm4 = F41 + F42 + F43 
# F41 = F41 / norm4
# F42 = F42 / norm4
# F43 = F43 / norm4


# Constructing Gebhart factors to estimate view factors of gray bodies (considering radiative rebounds)
# Note: matrix is constructed in blocks
Id44 = np.eye(4) # Identity matrix
Zero44 = np.zeros((4, 4)) # Zero matrix
# sparse matrix as we do not want 1 & 3 to see each other
Sparse13 = np.array([
     [1, 0, 0, 0],
     [0, 1, 0, 0],
     [0, 0, 0, 0],
     [0, 0, 0, 1]
     ])

Sparse31 = np.array([
     [0, 0, 0, 0],
     [0, 1, 0, 0],
     [0, 0, 1, 0],
     [0, 0, 0, 1]
     ])

lambda1 = -(1 - eps2) * F12
lambda2 = -(1 - eps1) * F21
lambda3 = -(1 - eps3) * F23
lambda4 = -(1 - eps2) * F32
lambda5 = -(1 - eps1) * F41
lambda6 = -(1 - eps2) * F42
lambda7 = -(1 - eps3) * F43

# Constructing the block matrix
Mat_A = np.block([
    [Id44, lambda1 * Sparse13, Zero44, Zero44],
    [lambda2 * Id44, Id44, lambda3 * Id44, Zero44],
    [Zero44, lambda4 * Sparse31, Id44, Zero44],
    [lambda5 * Id44, lambda6 * Id44, lambda7 * Id44, Id44]
])

# Matrix of known factors, simplified based on previous assumptions: 
# 1) Cryostat is a black body, 
# 2) panels 1 and 3 do not see each other even with GBVF, 
# 3) an object cannot see itself with VF.
Mat_EF = np.array([
    [0],
    [eps2 * F12],
    [0],
    [eps4 * F14],
    [eps1 * F21],
    [0],
    [eps3 * F23],
    [eps4 * F24],
    [0],
    [eps2 * F32],
    [0],
    [eps4 * F34],
    [eps1 * F41],
    [eps2 * F42],
    [eps3f * F43],
    [0]
])

# Solving the linear system Mat_A * BB = Mat_EF
BB = np.linalg.solve(Mat_A, Mat_EF)
# print('Matrix BB\n', BB)

# Displaying individual coefficients of the Gebhart factors
B11, B12, B13, B14, B21, B22, B23, B24, B31, B32, B33, B34, B41, B42, B43, B44 = BB[0:16]

# Renormalizing to avoid energy losses
for i in range(1, 4):
    normBi = sum(BB[(i-1)*4:(i)*4])
    for j in range(4):
        BB[(i-1)*4 + j] /= normBi

# Display of GBVF
for i in range(4):
    for j in range(4):
        print(f'B{i+1}{j+1} =', BB[i*4 + j])


# Extracting values of GR12, GR23, GR24, GR34, necessary radiative conductances for the problem
GR12 = A1 * eps1 * float(B12)
GR23 = A2 * eps3 * float(B23)
GR24 = A2 * eps2 * float(B24)
GR34 = A3 * eps3f * float(B34)
GR21 = A2 * eps2 * float(B21)
GR32 = A3 * eps3 * float(B32)
GR14 = A1 * eps1 * float(B14)

# Now seeking to find the temperature values of the panels
# System of 2 equations to solve: S * T = Z, derived from the node equation and the boundary temperature conditions obtained
Z = np.array([
    [0],
    [0],
    [T1**4],
    [Ts**4]
])

S = np.array([
    [GR21 * sigma,   -GR21 * sigma - GR23 * sigma - GR24 * sigma,     GR23 * sigma,   GR24 * sigma],
    [0,   GR32 * sigma,   -GR32 * sigma - GR34 * sigma,    GR34 * sigma],
    [1, 0, 0, 0],
    [0, 0, 0, 1]
])

T = np.linalg.solve(S, Z)
# print(T)
T = np.sqrt(np.sqrt(T)) # Convert from power to temperature
T1 = float(T[0])
T2 = float(T[1])
T3 = float(T[2])
T4 = float(T[3])
print('\n-----\nTemperatures of V Groove panels \nPanel 1:', T1, "K\nPanel 2:", T2, "K\nPanel 3:", T3, "K\nCryostat/Space:", T4, "K")

# Calculating radiative heat transfers
Q12_rad = sigma * eps1 * A1 * B12 * (T1**4 - T2**4)
Q14_rad = sigma * eps1 * A1 * B14 * (T1**4 - T4**4)
R_VG1S = Q12_rad / Q14_rad # Ratio of emission for panel 1
print('-----\nEmission ratio VG1', float(R_VG1S))
Q23_rad = sigma * eps2 * A2 * B23 * (T2**4 - T3**4)
Q24_rad = sigma * eps2 * A2 * B24 * (T2**4 - T4**4)
R_VG2S = Q23_rad / Q24_rad # Ratio of emission for panel 2
print('Emission ratio VG2', float(R_VG2S))
Q34_rad = sigma * eps3f * A3 * B34 * (T3**4 - T4**4)
P_cryo = int((Q14_rad + Q24_rad + Q34_rad) * 1000) # Power dissipated by the cryostat
print('-----\nPower dissipated by cryostat:', float(P_cryo), 'mW')

