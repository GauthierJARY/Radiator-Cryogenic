# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 15:40:15 2024

@author: jaryg
"""


import matplotlib.pyplot as plt
import matplotlib.patches as patches

def plot_summary_model_gauthier(T00, T0, T1, T2, T3, T4, Q_cond_01, Q_cond_04, Q_rad_001, Q_rad_004, Q_h0,
                                 Q_cond_12, Q_rad_12, Q_rad_14, Q_cond_23, Q_rad_23, Q_rad_24, Q_h3, Q_rad_34):

    # Calcul des températures et flux
    temp1 = int(T1)  # Température du bloc 1
    temp2 = int(T2)  # Température du bloc 2
    temp3 = int(T3)  # Température du bloc 3
    temp4 = T4  # Température du bloc 4
    temp5 = int(T0)
    temp6 = int(T00)
    
    flux_conduction_base_1 = int(float(Q_cond_01) * 1000)  
    flux_conduction_base_cryo = int(float(Q_cond_04) * 1000)  
    flux_radiation_MLI_1 = int(float(Q_rad_001) * 1000) 
    flux_radiation_MLI_cryo = int(float(Q_rad_004) * 1000) 
    flux_heater_base = int(float(Q_h0) * 1000) 
    
    flux_conduction_1_2 = int(float(Q_cond_12) * 1000)  
    flux_radiation_1_2 = int(float(Q_rad_12) * 1000)  
    flux_radiation_1_4 = int(float(Q_rad_14) * 1000)  
    
    flux_conduction_2_3 = int(float(Q_cond_23) * 1000)  
    flux_radiation_2_3 = int(float(Q_rad_23) * 10000) / 10  
    flux_radiation_2_4 = int(float(Q_rad_24) * 1000) 
      
    flux_heater_3 = int(float(Q_h3) * 1000) 
    flux_radiation_3_4 = int(float(Q_rad_34) * 1000)  
    
    # Création de la figure et des axes
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Création des blocs
    rect1 = patches.Rectangle((0, 3.15), 1, 0.5, facecolor='lightblue', edgecolor='black')
    rect2 = patches.Rectangle((0, 1.78), 1, 0.5, facecolor='lightgreen', edgecolor='black')
    rect3 = patches.Rectangle((0, 0.4), 1, 0.5, facecolor='lightcoral', edgecolor='black')
    rect4 = patches.Rectangle((0, -0.75), 1, 0.5, facecolor='lightyellow', edgecolor='red')  # Rectangle en dessous
    rect5 = patches.Rectangle((1.75, -1), 0.5, 5, facecolor='lightgray', edgecolor='black')  # Rectangle sur le côté
    
    # Ajout des blocs à l'axe
    ax.add_patch(rect1)
    ax.add_patch(rect2)
    ax.add_patch(rect3)
    ax.add_patch(rect4)
    ax.add_patch(rect5)
    
    # Ajout des noms et températures
    ax.text(0.5, 3.4, f"VG3\n {temp3}K", ha='center', va='center')
    ax.text(0.5, 2, f"VG2\n {temp2}K", ha='center', va='center')
    ax.text(0.5, 0.65, f"VG1\n {temp1}K", ha='center', va='center')
    ax.text(0.5, -0.5, f"Base\n {temp5}K & MLI {temp6}K", ha='center', va='center')  # Ajuster la position du texte
    ax.text(2, 1.5, f"Cryostat\n {temp4}K", ha='center', va='center')  # Ajuster la position du texte pour le rectangle sur le côté
    
    # Ajout des flèches de flux
    # Flux 2&3
    ax.annotate("", xy=(0.7, 2.3), xytext=(0.7, 3.2),
                arrowprops=dict(arrowstyle='<-', color='blue', lw=3))
    ax.text(0.85, 2.60, f"{flux_conduction_2_3}", ha='center', color='blue')
    
    ax.annotate("", xy=(0.3, 2.3), xytext=(0.3, 3.2),
                arrowprops=dict(arrowstyle='<-', color='orange', lw=3))
    ax.text(0.15, 2.6, f"{flux_radiation_2_3}", ha='center', color='orange')
    
    ax.annotate("", xy=(-0.5, 3.5), xytext=(0, 3.5),
                arrowprops=dict(arrowstyle='<-', color='red', lw=3))
    ax.text(-0.3, 3.3, f"{flux_heater_3 }", ha='center', color='red')
    
    ax.annotate("", xy=(1, 3.5), xytext=(1.75, 2.3),
                arrowprops=dict(arrowstyle='<-', color='orange', lw=3))
    ax.text(1.15, 3.5, f"{flux_radiation_3_4}", ha='center', color='orange')
    
    # Flux 1&2
    ax.annotate("", xy=(0.3, 0.9), xytext=(0.3, 1.8),
                arrowprops=dict(arrowstyle='<-', color='orange', lw=3))
    ax.text(0.15, 1.1, f"{flux_radiation_1_2}", ha='center', color='orange')
    
    ax.annotate("", xy=(0.7, 0.9), xytext=(0.7, 1.8),
                arrowprops=dict(arrowstyle='<-', color='blue', lw=3))
    ax.text(0.85, 1.1, f"{flux_conduction_1_2}", ha='center', color='blue')
    
    ax.annotate("", xy=(1, 2), xytext=(1.75, 1.8),
                arrowprops=dict(arrowstyle='<-', color='orange', lw=3))
    ax.text(1.15, 2.1, f"{flux_radiation_2_4}", ha='center', color='orange')
    
    ax.annotate("", xy=(1, 0.5), xytext=(1.75, 1.6),
                arrowprops=dict(arrowstyle='<-', color='orange', lw=3))
    ax.text(1.2, 0.5, f"{flux_radiation_1_4}", ha='center', color='orange')
    
    # Flux 0&1
    ax.annotate("", xy=(-0.5, -0.6), xytext=(0, -0.6),
                arrowprops=dict(arrowstyle='<-', color='red', lw=3))
    ax.text(-0.3, -0.8, f"{flux_heater_base}", ha='center', color='red')
    
    ax.annotate("", xy=(0.3, -0.25), xytext=(0.3, 0.45),
                arrowprops=dict(arrowstyle='<-', color='blue', lw=3))
    ax.text(0.15, -0.1, f"{flux_conduction_base_1 }", ha='center', color='blue')
    
    ax.annotate("", xy=(0.7, -0.25), xytext=(0.7, 0.45),
                arrowprops=dict(arrowstyle='<-', color='orange', lw=3))
    ax.text(0.85, -0.1, f"{flux_radiation_MLI_1 }", ha='center', color='orange')
    
    ax.annotate("", xy=(1, -0.5), xytext=(1.75, 1.0),
                arrowprops=dict(arrowstyle='<-', color='orange', lw=3))
    ax.text(1.10, -0.05, f"{flux_radiation_MLI_cryo }", ha='center', color='orange')
    
    ax.annotate("", xy=(1, -0.6), xytext=(1.75, -0.6),
                arrowprops=dict(arrowstyle='<-', color='blue', lw=3))
    ax.text(1.15, -0.8, f"{flux_conduction_base_cryo }", ha='center', color='blue')
    
    
    
    
    
    # Configuration de l'axe
    ax.set_xlim(-0.5, 3)
    ax.set_ylim(-1.5, 4)  # Ajuste les limites pour inclure le nouveau rectangle
    ax.axis('off')  # Masquer les axes
    
    # Affichage
    plt.title("Heat transfer with heaters load, on simplified model of 3 Vgroove [mW]")
    plt.show()

print("done !")
    