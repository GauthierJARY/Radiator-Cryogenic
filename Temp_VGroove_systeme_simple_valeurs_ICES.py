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

# Constantes et paramètres

sigma = 5.8e-8 # W/m2, Constante Stefan Boltzman
eps1 = 0.023 # emissivté du panneau  V Groove 1
eps2 = 0.023 # emissivté du panneau  V Groove 2
eps3 = 0.023 # emissivté du panneau  V Groove 3 côté intérieur
eps3f = 0.8 # emissivté du panneau  V Groove 3 côté extérieur
eps4 = 1 # emissivté du cryostat / assimilé à l'espace
theta = 6 # degrés
theta = theta/360 * 2 * m.pi # radians
theta1 = theta
theta2 = theta
theta3 = theta
A1 = 0.0268973 # m2, surface du panneau V Groove 1
A2 = 0.0255773 *2 # m2, surface du panneau V Groove 1
A3 = 0.0244109 *2# m2, surface du panneau V Groove 1

T1 = 245 # K, température du panneau V Groove 1
Ts = 0 # K, température du du cryostat / espace

# Facteur de vue / Form factor
# NB: c'est un raisonnement purement géométrique, sans rebonds de rayon

# vue du panneau 1
F12 = 1 - m.sin(theta1/2)
F14 = m.sin(theta1/2)
# vue du panneau 2
F21 = A1/A2 * F12
F24 = 2 * m.sin(theta2/2)
F23 = 1 - F24 - F21
# vue du panneau 3
F32 = F23*A2/A3 # NB: calculé par la formule de réflexitivité, mais si on mets une valeur faible ici et très élévée proche de 1 pour F34, on a un meilleur rapport des puissances transmises et dissipées
F34 = 1 - F32 
# vue du cryostat / espace
F41 = 0.0 # on néglige, d'après la formule de réciprocité car on a Space avec surface quasi infinie, très élevée en tout cas 
F42 = 0.0
F43 = 0.0
F3o4 = 1
# vue intérieure et extérieur du panneau 3
F3i2 = 1 - m.sin(theta/2)
F23i = 1 - m.sin(theta/2)

# NB: Si les facteurs de vue ne sont pas bons ni normalisés, il faut adapter cela. On normalise
norm1 = F12 + F14
F12 = F12/norm1
F14 = F14/norm1
# norm2 = F21 + F23 + F24
# F21 = F21/norm2
# F23 = F23/norm2
# F24 = F24/norm2
norm3 = F32 + F34
F32 = F32/norm3
F34 = F34/norm3
norm4 = F41 + F42 + F43 
# F41 = F41/norm4
# F42 = F42/norm4
# F43 = F43/norm4


# Construction des facteurs de Gebhart pour estimer les facteur de vue des corps gris (on prend en compte les rebonds radiatifs)
# NB: on construit la matrice par blocs
Id44 = np.eye(4)
Zero44 = np.zeros((4,4))
Sparse13 = np.array([
     [1,0,0,0],
     [0,1,0,0],
     [0,0,0,0],
     [0,0,0,1]
     ])

Sparse31 = np.array([
     [0,0,0,0],
     [0,1,0,0],
     [0,0,1,0],
     [0,0,0,1]
     ])

lambda1 = -(1-eps2)*F12
lambda2 = -(1-eps1)*F21
lambda3 = -(1-eps3)*F23
lambda4 = -(1-eps2)*F32
lambda5 = -(1-eps1)*F41
lambda6 = -(1-eps2)*F42
lambda7 = -(1-eps3)*F43

# Construction de la matrice par blocs
Mat_A = np.block([
    [Id44, lambda1 * Sparse13, Zero44, Zero44],
    [lambda2 * Id44, Id44, lambda3 * Id44, Zero44],
    [Zero44, lambda4 * Sparse31, Id44, Zero44],
    [lambda5 * Id44, lambda6 * Id44, lambda7 * Id44, Id44]
])

# matrice des facteurs connus, et simplifiés en fonction des hypothèses préalable: 
# 1) cryostat corps noir, 
# 2) 1 et 3 ne se voient pas, 
# 3) un objet ne peut pas se voir lui même.
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

# Résolution du système linéaire Mat_A * BB = Mat_EF
BB = np.linalg.solve(Mat_A, Mat_EF)
# print('Matrice BB\n',BB)

# Affichage des coefficients individuels des facteurs de Gebhart
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

# On renormalise pour eviter les pertes d'énergie
normB1 = B11 + B12 + B13 + B14
B11 = B11 / normB1
B12 = B12 / normB1
B13 = B13 / normB1
B14 = B14 / normB1
normB2 = B21 + B22 + B23 + B24
B21 = B21 / normB2
B22 = B22 / normB2
B23 = B23 / normB2
B24 = B24 / normB2
normB3 = B31 + B32 + B33 + B34
B31 = B31 / normB3
B32 = B32 / normB3
B33 = B33 / normB3
B34 = B34 / normB3

# print('B11 =', B11)
# print('B12 =', B12)
# print('B13 =', B13)
# print('B14 =', B14)
# print('B21 =', B21)
# print('B22 =', B22)
# print('B23 =', B23)
# print('B24 =', B24)
# print('B31 =', B31)
# print('B32 =', B32)
# print('B33 =', B33)
# print('B34 =', B34)
# print('B41 =', B41)
# print('B42 =', B42)
# print('B43 =', B43)
# print('B44 =', B44)
# Extraction des valeurs de GR12, GR23, GR24, GR34, conductances radiatives nécessaires de notre problèmes
GR12 = A1 * eps1 * float(B12)
GR23 = A2 * eps3 * float(B23)
GR24 = A2 * eps2 * float(B24)
GR34 = A3 * eps3f * float(B34)
GR21 = A2 * eps2 * float(B21)
GR32 = A3 * eps3 * float(B32)
GR14 = A1 * eps1 * float(B14)


# On cherche maintenant à trouver les valeurs des températures des panneaux
# Système de 2 équations à résoudre : S * T = Z, issu de l'équation aux noeuds, et des 2 conditions limites en température que l'on a obtenu
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
# print(T)
T = np.sqrt(np.sqrt(T))
T1 = float(T[0])
T2 = float(T[1])
T3 = float(T[2])
T4 = float(T[3])
print('\n-----\nTemperatures des panneaux V Groove \nPanneau 1:',T1, "K\nPanneau 2:",T2, "K\nPanneau 3:",T3, "K\nCryostat/ Espace:",T4, "K")

Q12_rad = sigma * eps1 * A1 * B12 * (T1**4 - T2**4)
Q14_rad = sigma * eps1 * A1 * B14 * (T1**4 - T4**4)
R_VG1S = Q12_rad/Q14_rad
print('-----\nRapport emission VG1', float(R_VG1S))
Q23_rad = sigma * eps2 * A2 * B23 * (T2**4 - T3**4)
Q24_rad = sigma * eps2 * A2 * B24 * (T2**4 - T4**4)
R_VG2S = Q23_rad/Q24_rad
print('Rapport emission VG2', float(R_VG2S))
Q34_rad = sigma * eps3f * A3 * B34 * (T3**4 - T4**4)
P_cryo =int((Q14_rad + Q24_rad + Q34_rad)*1000)
print('-----\nPuissance dissipée par cryostat:', float(P_cryo),'mW')


#### IGNORE as we have correct results above
#### Autre methode de calcul de gebhart one by one 
# Résolution du système linéaire S * T = Z
Mat_C = np.array([
    [1, -(1-eps2)*F12, 0, 0], 
    [-(1-eps1)*F21, 1, -(2-eps3)*F23, 0], 
    [0, -(1-eps2)*F32, 1, 0], 
    [-(1-eps1)*F41, -(1-eps2)*F42, -(1-eps3)*F43, 1]
    ])

Y = np.array([[0], [eps1*F21], [0], [eps1*F41]])

# Inversion de la matrice Mat_C
Mat_C_inv = np.linalg.inv(Mat_C)

# Résolution de l'équation Mat_C * X = Y
X = np.dot(Mat_C_inv, Y)

# Affichage du résultat
# print('Solution X =\n', X)
