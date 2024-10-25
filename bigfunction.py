# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 14:08:29 2024

@author: jaryg
"""

from scipy.optimize import fsolve
import math as m
import numpy as np
import matplotlib.pyplot as plt
from GreyBodyViewFactor_3VG_analytical import calculate_Gebharts

# copy paste this in console
# cd C:\Users\jaryg\Desktop\Tests Python



# Variables
Q_heater_base = 0.85 # W
T0_values =[]
T00_values =[]
T1_values =[]
T2_values =[]
T3_values =[]
T4_values =[]
Q_heater_base_values = np.linspace(5,12,10)/10

for Q_heater_base in Q_heater_base_values:
    Q_heater_VG3 = 50e-3 #W
    L_cryo = 0.710/2 # length rod of fixation to cryostat
    d_steel = 5e-3 # diameter rod of fixation to cryostat
    # Constants
    sigma = 5.8e-8 # W/m2, Constante Stefan Boltzman
    k_nylon = 0.25
    k_steel = 12
    # Parameters
    T4 = 4 # K, temperature cryostat
    L_nylon_screw = 8e-3 #m
    d_ext = 0.01
    d_int = 0.004
    S_nylon_screw = m.pi/4*(d_ext**2-d_int**2)
        # VGroove
    eps1 = 0.05 # emissivté du panneau  V Groove 1
    eps2 = 0.05 # emissivté du panneau  V Groove 2
    eps3 = 0.05 # emissivté du panneau  V Groove 3 côté intérieur
    eps3f = 0.7 # emissivté du panneau  V Groove 3 côté extérieur
    eps4 = 1 # emissivté du cryostat / assimilé à l'espace
    theta = 5 # degrés
    theta = theta/360 * 2 * m.pi # radians
    theta1 = theta
    theta2 = theta
    theta3 = theta
    A_VG1 = 0.0268973 # m2, surface du panneau V Groove 1
    A_VG2 = 0.0255773 # m2, surface du panneau V Groove 1
    A_VG3 = 0.0244109 # m2, surface du panneau V Groove 1
        # Model of MLI
    eps_MLI = 0.05
    A_MLI =  0.115 #m2
    S_steel= m.pi/4 * d_steel**2
    
    k_eff_MLI = 0.05 # cf Lornenzo mail
    S_MLI = 0.115 # m2
    eps_eff_MLI = 0.028 # cf Lornenzo mail
    eps_MLI = 0.05

    
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
    
    B12 = B12/2 # il voit désormais la base et pas l'espace
    B100 = B12
    B001 = B100 * A_VG1/S_MLI # on utilise la formule de réciprocité
    B004 = 1-B001
    
    GR12 = A_VG1 * eps1 * float(B12)
    GR23 = A_VG2 * eps2 * float(B23)
    GR24 = A_VG2 * eps2 * float(B24)
    GR34 = A_VG3 * eps3f * float(B34)
    GR21 = A_VG2 * eps2 * float(B21)
    GR32 = A_VG3 * eps3 * float(B32)
    GR14 = A_VG1 * eps1 * float(B14)
    GR12 = max(GR12,GR21)
    GR21 = max(GR12,GR21)
    GR23 = max(GR23,GR32)
    GR32 = max(GR23,GR32)
    

    GL12 = k_nylon/L_nylon_screw*S_nylon_screw
    GL23 = k_nylon/L_nylon_screw*S_nylon_screw
    
    GL04 = k_steel*S_steel/L_cryo
    GL01 = k_nylon/L_nylon_screw*S_nylon_screw
    GL000 = k_eff_MLI * S_MLI 
    
    GR000 = eps_eff_MLI * S_MLI
    GR004 = eps_MLI * S_MLI * B004
    GR001 = eps_MLI * S_MLI * B001
    
    def prediction_temperatures_from_heater_gauthier(vars):
        T0x, T00x, T1x, T2x, T3x = vars
        eq1 = - Q_heater_base - GL01*(T1x-T0x) - GL000*(T00x-T0x) - GL04*(T4-T0x) + sigma * GR000 *(T0x**4-T00x**4)
        eq2 =  -(-GL000*(T00x-T0x)) - sigma * GR000 *(T0x**4-T00x**4) + sigma * GR004 *(T00x**4-T4**4) + sigma * GR001*(T00x**4-T1x**4) 
        eq3 = -(- GL01 * (T1x - T0x)) - GR001 * sigma * (T00x**4 - T1x*4) + GR14 * sigma * (T1x**4 - T4**4) - GL12 * (T2x - T1x) + GR12 * sigma * (T1x**4 - T2x**4)
        eq4 = -(- GL12 * (T2x - T1x)) - GR12 * sigma * (T1x**4 - T2x**4) + GR23 * sigma * (T2x**4 - T3x**4) - GL23 * (T3x - T2x) + GR24 * sigma * (T2x**4- T4**4)
        eq5 = - GR32 * sigma * (T2x**4- T3x**4) + GR34 * sigma * (T3x**4- T4**4) -(- GL23 * (T3x - T2x)) - Q_heater_VG3
        
        return [eq1, eq2, eq3, eq4, eq5]
    
    T0, T00, T1, T2, T3 =  fsolve(prediction_temperatures_from_heater_gauthier, (320, 250, 200, 120, 100))
    print(f'\n-----\nTemperatures of V Groove \nPannel 1: {T1:.1f} K\nPannel 2: {T2:.1f} K\nPannel 3: {T3:.1f} K\nCryostat/ Space: {T4:.1f} K\nMLI: {T00:.1f} K\nBase: {T0:.1f} K')
    
    Q_cond_01  = - GL01 * (T1 - T0)
    Q_cond_04 = - GL04*(T4-T0)
    Q_cond_000 = - GL000*(T00-T0)
    Q_rad_001 = sigma * GR001*(T00**4-T1**4)
    Q_rad_004 = sigma * GR004 *(T00**4-T4**4)
    Q_rad_000 = sigma * GR000*(T0**4-T00**4)
    Q_h0 = Q_heater_base
    Q_cond_12 = - GL12 * (T2 - T1)
    Q_rad_12 = GR21 * sigma * (T1**4 - T2**4)
    Q_rad_14 = GR14 * sigma * (T1**4 - T4**4)
    Q_cond_23 = - GL23 * (T3 - T2)
    Q_rad_23 = GR23 * sigma * (T2**4 - T3**4)
    Q_rad_24 = GR24 * sigma * (T2**4- T4**4)
    Q_h3 = Q_heater_VG3
    Q_rad_34= GR34 * sigma * (T3**4- T4**4)
    
    
    from plot_functions import plot_summary_model_gauthier
    
    plot_summary_model_gauthier(T00, T0, T1, T2, T3, T4, Q_cond_01, Q_cond_04, Q_rad_001, Q_rad_004, 
                                  Q_h0, Q_cond_12, Q_rad_12, Q_rad_14, Q_cond_23, Q_rad_23, 
                                  Q_rad_24, Q_h3, Q_rad_34)
    T0_values.append(T0)
    T00_values.append(T00)
    T1_values.append(T1)
    T2_values.append(T2)
    T3_values.append(T3)

# Tracer les courbes
plt.figure(figsize=(10, 6))
plt.plot(Q_heater_base_values, T0_values, label='T0')
plt.plot(Q_heater_base_values, T00_values, label='T00')
plt.plot(Q_heater_base_values, T1_values, label='T1')
plt.plot(Q_heater_base_values, T2_values, label='T2')
plt.plot(Q_heater_base_values, T3_values, label='T3')

# Ajout de labels et titre
plt.title("temperatures vs load from base heater")
plt.xlabel("Q_heater_base input")
plt.ylabel("Temperature (K)")
plt.legend()
plt.grid()


