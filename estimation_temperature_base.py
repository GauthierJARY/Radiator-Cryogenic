# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 11:26:32 2024

@author: jaryg

Objectif: estimer la charge en puissance d'un réchauffeur installé sur la base de la maquette. 
Celle ci est attachée au socle du croystat, doit avoir une certaine température et le cryostat dispose 
d'une certaine puissance froide maximum'
"""

import numpy as np
import math as m
import matplotlib.pyplot as plt

T_cryo = 4 #K
T_base = 300 #K
T_VG1 = 183 #K
sigma = 5.67e-8
eps_MLI = 0.05
A_base =  0.076 #m2
k_steel = 12
k_nylon = 0.29
L_cryo = 0.710/2
L_nylon = 0.0075
d_vis = 4e-3
d_steel = 5e-3
nb_vis = 3
S_nylon = nb_vis * m.pi/4 * d_vis**2
S_steel= m.pi/4 * d_steel**2
Gebhart_base_VG1 = 0.3
GR_base_VG1 = eps_MLI * A_base * Gebhart_base_VG1
Gebhart_base_cryo = 0.7
GR_base_cryo = eps_MLI * A_base * Gebhart_base_cryo


## We do a loop on the temperature to see if a lower temperature of the base would reduce the power needed and its sensibility
# We do equilibrium of heat exchange for several value of T_base
# In addition, we also seek to investigate impact of the rod diameter if it is made out of steel 
list_temperature_base =np.linspace(230,300,10)
list_diameter_rod = np.linspace(5,15,5)*1e-3
for d_steel in list_diameter_rod:   
    S_steel= m.pi/4 * d_steel**2
    
    list_power =[]
    for T_base in list_temperature_base: 
        Q_rad_cryo = sigma*GR_base_cryo*(T_cryo**4 - T_base**4)
        Q_rad_VG1 = sigma*GR_base_VG1*(T_VG1**4 - T_base**4)
        Q_cond_VG1 = k_nylon*S_nylon/L_nylon*(T_VG1 - T_base)
        Q_cond_cryo = k_steel*S_steel/L_cryo*(T_cryo - T_base)
        Q_heater =  Q_rad_cryo +Q_cond_cryo + Q_rad_VG1 + Q_cond_VG1 
        
        print('#########################################\n##  Power of Heater needed', int(100000*Q_heater)/100, 'mW  ##\n#########################################')
        print('Radiative Flux to VG1', int(100000*Q_rad_VG1)/100, 'mW')
        print('Radiative Flux to cryo', int(100000*Q_rad_cryo)/100, 'mW')
        print('Conductive Flux to VG1', int(100000*Q_cond_VG1)/100, 'mW')
        print('Conductive Flux to cryo', int(100000*Q_cond_cryo)/100, 'mW')
        print('Total Flux to cryo', int(100000*(Q_cond_cryo + Q_rad_cryo))/100, 'mW')
        list_power.append(-Q_heater*1000)
    plt.plot(list_temperature_base,list_power, label=f'Rod diameter {1000*d_steel} mm')

plt.title('Power of the heater to reach base temperature')
plt.xlabel('Base temperature [K]')
plt.ylabel('Power [mW]')
plt.yticks(np.arange(400, 4001, 250))
plt.xticks(np.arange(230, 301, 10)) 
plt.legend()
plt.grid()
plt.show()

## We do a loop on the temperature to see if a lower temperature of the base would reduce the power needed and its sensibility
# We do equilibrium of heat exchange for several value of T_base
# In addition, we also seek to investigate impact of the temperature of VG1

list_temperature_base =np.linspace(230,300,10)
list_temperature_VG1 =np.linspace(150,250,10)
d_steel = 5e-3
S_steel= m.pi/4 * d_steel**2
for T_VG1 in list_temperature_VG1:      
    list_power =[]
    for T_base in list_temperature_base: 
        Q_rad_cryo = sigma*GR_base_cryo*(T_cryo**4 - T_base**4)
        Q_rad_VG1 = sigma*GR_base_VG1*(T_VG1**4 - T_base**4)
        Q_cond_VG1 = k_nylon*S_nylon/L_nylon*(T_VG1 - T_base)
        Q_cond_cryo = k_steel*S_steel/L_cryo*(T_cryo - T_base)
        Q_heater =  Q_rad_cryo +Q_cond_cryo + Q_rad_VG1 + Q_cond_VG1 
        
        print('#########################################\n##  Power of Heater needed', int(100000*Q_heater)/100, 'mW  ##\n#########################################')
        print('Radiative Flux to VG1', int(100000*Q_rad_VG1)/100, 'mW')
        print('Radiative Flux to cryo', int(100000*Q_rad_cryo)/100, 'mW')
        print('Conductive Flux to VG1', int(100000*Q_cond_VG1)/100, 'mW')
        print('Conductive Flux to cryo', int(100000*Q_cond_cryo)/100, 'mW')
        print('Total Flux to cryo', int(100000*(Q_cond_cryo + Q_rad_cryo))/100, 'mW')
        list_power.append(-Q_heater*1000)
    plt.plot(list_temperature_base,list_power, label=f'Temp VG1 {int(T_VG1)} K')

plt.title('Power of the heater to reach base temperature')
plt.xlabel('Base temperature [K]')
plt.ylabel('Power [mW]')
plt.yticks(np.arange(400, 2001, 250))
plt.xticks(np.arange(230, 301, 10)) 
plt.legend()
plt.grid()
plt.show()

## We do a loop on the temperature to see if a lower temperature of the base would reduce the power needed and its sensibility
# We do equilibrium of heat exchange for several value of T_base
# In addition, we also seek to investigate impact of the emissivity of the MLI, and seems to be the most determinant criteria

list_temperature_base =np.linspace(230,300,10)
list_emissivity_MLI = np.linspace(1,5,5)/100
for eps_MLI in list_emissivity_MLI:  
    GR_base_VG1 = eps_MLI * A_base * Gebhart_base_VG1
    GR_base_cryo = eps_MLI * A_base * Gebhart_base_cryo    
    list_power =[]
    for T_base in list_temperature_base: 
        Q_rad_cryo = sigma*GR_base_cryo*(T_cryo**4 - T_base**4)
        Q_rad_VG1 = sigma*GR_base_VG1*(T_VG1**4 - T_base**4)
        Q_cond_VG1 = k_nylon*S_nylon/L_nylon*(T_VG1 - T_base)
        Q_cond_cryo = k_steel*S_steel/L_cryo*(T_cryo - T_base)
        Q_heater =  Q_rad_cryo +Q_cond_cryo + Q_rad_VG1 + Q_cond_VG1 
        
        print('#########################################\n##  Power of Heater needed', int(100000*Q_heater)/100, 'mW  ##\n#########################################')
        print('Radiative Flux to VG1', int(100000*Q_rad_VG1)/100, 'mW')
        print('Radiative Flux to cryo', int(100000*Q_rad_cryo)/100, 'mW')
        print('Conductive Flux to VG1', int(100000*Q_cond_VG1)/100, 'mW')
        print('Conductive Flux to cryo', int(100000*Q_cond_cryo)/100, 'mW')
        print('Total Flux to cryo', int(100000*(Q_cond_cryo + Q_rad_cryo))/100, 'mW')
        list_power.append(-Q_heater*1000)
    plt.plot(list_temperature_base,list_power, label=f'Emissivity MLI {eps_MLI}')

plt.title('Power of the heater to reach base temperature')
plt.xlabel('Base temperature [K]')
plt.ylabel('Power [mW]')
plt.yticks(np.arange(400, 2001, 250))
plt.xticks(np.arange(230, 301, 10)) 
plt.legend()
plt.grid()
plt.show()

# # quick conclusion: 
#     we can play on some paramters, 
#     -emissivity of MLI (we can do a good MLI)
#     -length of steel rod 
#     -conductivity of fixation rod
#     -diameter of fixation rod
#     -temperature of the base
# # CONCLUSION
# So the idea would now be to assess a design final on NX Siemens CAD, 
# put boundaries and try to assess the input power of the heater that 
# we would need with analytical method and then try and error for optimisation (script ?)
# THE MAIN THING TO WORK ON IS MLI efficiency
# The second main restriction will first be structural : the rod need to support 
# all model and instrumentation.

# Optimized one? 
print('\n\n ### OPTIMIZED###')
T_cryo = 4 #K
T_base = 290 #K
T_VG1 = 183 #K
sigma = 5.67e-8
eps_MLI = 0.01  # good MLI optimized
A_base =  0.076 #m2
k_steel = 12 # we put nylon instead, and bigger rod
k_nylon = 0.29
L_cryo = 0.710*2/3
L_nylon = 0.0075
d_vis = 4e-3
d_steel = 1e-2 # bigger rod since it is from nylon
nb_vis = 3
S_nylon = nb_vis * m.pi/4 * d_vis**2
S_steel= m.pi/4 * d_steel**2
Gebhart_base_VG1 = 0.3
GR_base_VG1 = eps_MLI * A_base * Gebhart_base_VG1
Gebhart_base_cryo = 0.7
GR_base_cryo = eps_MLI * A_base * Gebhart_base_cryo

Q_rad_cryo = sigma*GR_base_cryo*(T_cryo**4 - T_base**4)
Q_rad_VG1 = sigma*GR_base_VG1*(T_VG1**4 - T_base**4)
Q_cond_VG1 = k_nylon*S_nylon/L_nylon*(T_VG1 - T_base)
Q_cond_cryo = k_steel*S_steel/L_cryo*(T_cryo - T_base)
Q_heater =  Q_rad_cryo +Q_cond_cryo + Q_rad_VG1 + Q_cond_VG1 

print('#########################################\n##  Power of Heater needed', int(100000*Q_heater)/100, 'mW  ##\n#########################################')
print('Radiative Flux to VG1', int(100000*Q_rad_VG1)/100, 'mW')
print('Radiative Flux to cryo', int(100000*Q_rad_cryo)/100, 'mW')
print('Conductive Flux to VG1', int(100000*Q_cond_VG1)/100, 'mW')
print('Conductive Flux to cryo', int(100000*Q_cond_cryo)/100, 'mW')
print('Total Flux to cryo', int(100000*(Q_cond_cryo + Q_rad_cryo))/100, 'mW')

# Données pour le camembert
labels = [ 'Radiative Flux VG1', 'Radiative Flux Cryo', 'Conductive Flux VG1', 'Conductive Flux Cryo']
sizes = [
    -Q_rad_VG1,
    -Q_rad_cryo,
    -Q_cond_VG1,
    -Q_cond_cryo
]
# Normalisation des valeurs pour qu'elles somment à 100
total = sum(sizes)
sizes = [(size / total) * 100 for size in sizes]
# Création du diagramme en camembert
plt.figure(figsize=(8, 8))
wedges, texts, autotexts = plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140)
# Rapprochement des labels
# Personnalisation des étiquettes et pourcentages
for i, text in enumerate(texts):
    text.set_color(wedges[i].get_facecolor())  # Couleur des labels
    text.set_fontsize(12)                       # Taille des labels
    text.set_fontweight('bold')                  # Gras
    text.set_position((0.9 * text.get_position()[0], 1 * text.get_position()[1]))  # Décalage vers l'extérieur
for autotext in autotexts:
    autotext.set_fontsize(14)                   # Taille des pourcentages
    autotext.set_fontweight('bold')              # Gras
plt.setp(autotexts, color='white')             # Couleur des pourcentages en blanc
plt.title('Flux Thermique', fontsize=16, fontweight='bold')  # Titre
plt.axis('equal')  # Assure que le camembert est bien circulaire
plt.show()