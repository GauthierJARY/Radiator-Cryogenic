# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 16:42:47 2024

@author: jaryg
"""


import numpy as np
import math as m
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

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



# Paramètres constants
k_eff_MLI = 0.05      # Remplacez par votre valeur
S_MLI = A_base 
eps_eff_MLI = 0.028 # Remplacez par votre valeur
sigma = 5.67e-8  # Constante de Stefan-Boltzmann en W/(m² K⁴)

# Fonction à résoudre
def equation(T_MLI):
    return k_eff_MLI * S_MLI * (T_MLI - T_base) - (
        eps_eff_MLI * S_MLI * sigma * (T_VG1**4 - T_MLI**4) + 
        eps_eff_MLI * S_MLI * sigma * (T_cryo**4 - T_MLI**4)
    )
# Estimation initiale pour T_MLI
T_MLI_initial_guess = T_base  # Remplacez par une estimation raisonnable

# Résolution de l'équation
T_MLI = int(fsolve(equation, T_MLI_initial_guess))


# Optimized one? 
print('\n\n ### OPTIMIZED###')
T_cryo = 4 #K
T_base = 290 #K
T_VG1 = 183 #K
sigma = 5.67e-8
eps_MLI = 0.01  # good MLI optimized
A_base =  0.115 #m2
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
T_MLI = int(fsolve(equation, T_MLI_initial_guess))

Q_cond_MLI = k_eff_MLI * S_MLI * (T_MLI - T_base)
Q_cond_VG1 = k_nylon*S_nylon/L_nylon*(T_VG1 - T_base)
Q_cond_cryo = k_steel*S_steel/L_cryo*(T_cryo - T_base)
Q_heater =  Q_cond_MLI + Q_cond_cryo + Q_cond_VG1 

print('#########################################\n##  Power of Heater needed', int(100000*Q_heater)/100, 'mW  ##\n#########################################')
print('Conductive Flux to MLI', int(100000*Q_cond_MLI)/100, 'mW')
print('Conductive Flux to VG1', int(100000*Q_cond_VG1)/100, 'mW')
print('Conductive Flux to cryo', int(100000*Q_cond_cryo)/100, 'mW')
print('Total Flux to cryo', int(100000*(Q_cond_cryo))/100, 'mW')

# Données pour le camembert
labels = [ 'Conductive Flux MLI', 'Conductive Flux VG1', 'Conductive Flux Cryo']
sizes = [
    -Q_cond_MLI,
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