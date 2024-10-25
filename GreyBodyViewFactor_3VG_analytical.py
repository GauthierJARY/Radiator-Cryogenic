# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 16:39:48 2024

@author: jaryg
"""

import numpy as np
import math as m


def calculate_Gebharts(eps1, eps2, eps3,eps3f, eps4, theta, A_VG1, A_VG2, A_VG3):
    # We get the grey body view factor between VG from analytical
    # vue du panneau 1
    F12 = 1 - m.sin(theta/2)
    F14 = 1-F12
    # vue du panneau 2
    F21 = (1 - m.sin(theta/2))/2
    F23 = (1 - m.sin(theta/2))/2
    F24 = 1 - F21 - F23
    # vue du panneau 3
    F32 = F23*A_VG2/A_VG3 # NB: calculé par la formule de réflexitivité, mais si on mets une valeur faible ici et très élévée proche de 1 pour F34, on a un meilleur rapport des puissances transmises et dissipées
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
    # Extraction des coefficients
    B_values = BB.flatten()

    # Normalisation des coefficients
    B_normalized = []
    for i in range(4):
        norm = sum(B_values[i * 4:(i + 1) * 4])
        B_normalized.append(B_values[i * 4:(i + 1) * 4] / norm)
        
        
    # B_normalized[1] = B_normalized[1]/2
    return np.array(B_normalized).flatten()

print("done")

