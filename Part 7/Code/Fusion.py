# Regular imports
import numpy             as np
import matplotlib.pyplot as plt
import pickle as pkl 

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

#Unpickling launch
with open("Mission.pkl", 'rb') as file:
        mission = pkl.load(file)

"""
This code computes the answer to D.2)
"""

T_r = mission.system.star_temperature
r   = mission.system.star_radius * 1000
rho = (mission.system.star_mass * const.m_sun) / ((4/3)*np.pi*(r**3))

G   = const.G
m_h = const.m_p
k   = const.k_B

# We assume mean molecular mass to be 1
mu  = 1.75

T_c = T_r + (2/3)*np.pi*G*rho*(r**2)*(mu*m_h/k)

print(f"Core temperature : {T_c:.5e} K")

"""
This code computes the answer to D.3)
"""

# Computing reaction rate for the pp cycle
T_6   = T_c / (1e6)
e_0pp = 1.08 * 10**(-12)
X_H   = (0.745 * 1)/(0.745 * 1 + 0.253 * 4 + 0.002 * ((12 + 13 + 15)/3))
e_pp  = e_0pp * rho * (X_H**2) * (T_6**4)

# Computing reaction rate for the CNO
e_0CNO = 8.24 * 10**(-31)
X_CNO  = (0.002 * ((12 + 13 + 15)/3)) / (0.745 * 1 + 0.253 * 4 + 0.002 * ((12 + 13 + 15)/3)) 
e_CNO  = e_0CNO * X_H * X_CNO * rho* T_6**(20)

# Summing these rates together
e = e_pp + e_CNO
# Computing luminosity
L = e * (rho * (4/3) * np.pi * (0.2 * r)**3)
print(e)
print(f"Star luminosity  : {L/const.L_sun:.5e} W")

"""
Output :

Core temperature : 1.61498e+07 K
Star luminosity  : 5.26695e+23 W
"""