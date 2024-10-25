# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:51:44 2024

@author: jaryg
"""



import numpy as np
import math as m
import matplotlib.pyplot as plt
from GreyBodyViewFactor_3VG_analytical import calculate_Gebharts

# Constantes et paramètres
    
sigma = 5.8e-8 # W/m2, Constante Stefan Boltzman
eps1 = 0.05 # emissivté du panneau  V Groove 1
eps2 = 0.05 # emissivté du panneau  V Groove 2
eps3 = 0.05 # emissivté du panneau  V Groove 3 côté intérieur
eps3f = 0.7 # emissivté du panneau  V Groove 3 côté extérieur
eps4 = 1 # emissivté du cryostat / assimilé à l'espace

A_VG1 = 0.0268973 # m2, surface du panneau V Groove 1
A_VG2 = 0.0255773 # m2, surface du panneau V Groove 1
A_VG3 = 0.0244109 # m2, surface du panneau V Groove 1

T1 = 200 # K, température du panneau V Groove 1
Ts = 4 # K, température du du cryostat / espace
T1_values =[]
T2_values =[]
T3_values =[]
angle_values = np.linspace(0,150,25) 

for theta in angle_values:
    theta = theta/360 * 2 * m.pi # radians
    theta1 = theta
    theta2 = theta
    theta3 = theta      
        
    BB = calculate_Gebharts(eps1, eps2, eps3, eps3f, eps4, theta, A_VG1, A_VG2, A_VG3)
    B11 = BB[0]
    B12 = BB[1]
    B13 = BB[2]
    B14 = BB[3]
    B21 = BB[4]
    B22 = BB[5]
    B23 = BB[6]
    B24 = BB[7]
    B31 = BB[8]
    B32 = BB[9]
    B33 = BB[10]
    B34 = BB[11]
    B41 = BB[12]
    B42 = BB[13]
    B43 = BB[14]
    B44 = BB[15]
    
    GR12 = A_VG1 * eps1 * float(B12)
    GR23 = A_VG2 * eps2 * float(B23)
    GR24 = A_VG2 * eps2 * float(B24)
    GR34 = A_VG3 * eps3f * float(B34)
    GR21 = A_VG2 * eps2 * float(B21)
    GR32 = A_VG3 * eps3 * float(B32)
    GR14 = A_VG1 * eps1 * float(B14)
    Z = np.array([
        [0],
        [0],
        [T1**4],
        [Ts**4]
    ])
    S = np.array([
        [GR21 * sigma,   -GR21 * sigma -GR23 * sigma -GR24 * sigma,     GR23 * sigma,   GR24 * sigma],
        [0,   GR32 * sigma,   -GR32 * sigma - GR34 * sigma,    GR34 * sigma],
        [1,0,0,0],
        [0,0,0,1]
    ])

    T = np.linalg.solve(S, Z)
    T = np.sqrt(np.sqrt(T))
    T1 = float(T[0])
    T2 = float(T[1])
    T3 = float(T[2])
    T4 = float(T[3])
    T1_values.append(T1)
    T2_values.append(T2)
    T3_values.append(T3)

plt.figure(figsize=(10, 6))
plt.plot(angle_values, T1_values, label='T1')
plt.plot(angle_values, T2_values, label='T2')
plt.plot(angle_values, T3_values, label='T3')
plt.title("temperatures vs angle")
plt.xlabel("angle")
plt.ylabel("Temperature (K)")
plt.legend()
plt.grid()



sigma = 5.8e-8 # W/m2, Constante Stefan Boltzman
A_VG1 = 0.0268973 # m2, surface du panneau V Groove 1
A_VG2 = 0.0255773 # m2, surface du panneau V Groove 1
A_VG3 = 0.0244109 # m2, surface du panneau V Groove 1
theta = 6
theta = theta/360 * 2 * m.pi # radians
theta1 = theta
theta2 = theta
theta3 = theta    
T1 = 200 # K, température du panneau V Groove 1
Ts = 4 # K, température du du cryostat / espace
T1_values =[]
T2_values =[]
T3_values =[]
emissivity_values = np.linspace(1,10,10)/100
eps3f = 0.7 # emissivté du panneau  V Groove 3 côté extérieur
eps4 = 1 # emissivté du cryostat / assimilé à l'espace

for eps in emissivity_values:
    eps1 = eps 
    eps2 = eps
    eps3 = eps 
            
    BB = calculate_Gebharts(eps1, eps2, eps3, eps3f, eps4, theta, A_VG1, A_VG2, A_VG3)
    B11 = BB[0]
    B12 = BB[1]
    B13 = BB[2]
    B14 = BB[3]
    B21 = BB[4]
    B22 = BB[5]
    B23 = BB[6]
    B24 = BB[7]
    B31 = BB[8]
    B32 = BB[9]
    B33 = BB[10]
    B34 = BB[11]
    B41 = BB[12]
    B42 = BB[13]
    B43 = BB[14]
    B44 = BB[15]
    
    GR12 = A_VG1 * eps1 * float(B12)
    GR23 = A_VG2 * eps2 * float(B23)
    GR24 = A_VG2 * eps2 * float(B24)
    GR34 = A_VG3 * eps3f * float(B34)
    GR21 = A_VG2 * eps2 * float(B21)
    GR32 = A_VG3 * eps3 * float(B32)
    GR14 = A_VG1 * eps1 * float(B14)
    Z = np.array([
        [0],
        [0],
        [T1**4],
        [Ts**4]
    ])
    S = np.array([
        [GR21 * sigma,   -GR21 * sigma -GR23 * sigma -GR24 * sigma,     GR23 * sigma,   GR24 * sigma],
        [0,   GR32 * sigma,   -GR32 * sigma - GR34 * sigma,    GR34 * sigma],
        [1,0,0,0],
        [0,0,0,1]
    ])
    
    T = np.linalg.solve(S, Z)
    T = np.sqrt(np.sqrt(T))
    T1 = float(T[0])
    T2 = float(T[1])
    T3 = float(T[2])
    T4 = float(T[3])
    T1_values.append(T1)
    T2_values.append(T2)
    T3_values.append(T3)

plt.figure(figsize=(10, 6))
plt.plot(emissivity_values, T1_values, label='T1')
plt.plot(emissivity_values, T2_values, label='T2')
plt.plot(emissivity_values, T3_values, label='T3')
plt.title("temperatures vs emissivity")
plt.xlabel("emissivtiy")
plt.ylabel("Temperature (K)")
plt.legend()
plt.grid()