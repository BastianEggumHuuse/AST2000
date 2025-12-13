import numpy as np
import matplotlib.pyplot as plt

import ast2000tools.utils as utils
from   ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as const
from   ast2000tools.relativity import RelativityExperiments

from RelativeTrilateration import TrilaterationAlgorithm

LorentzFactor   = lambda v : 1/(1-v**2)

KgToMeters      = lambda kg : kg * (const.G / (const.c**2))

SecondsToMeters = lambda s : s*const.c


print("1)\n")

r_p = 2304.5943016 # in km
m_p = 2.754441653063e23

r_vec_1 = np.array([-3126.903,-4994.661])
h = np.linalg.norm(r_vec_1) - r_p

print(f"Height : {h} km\n")

print("2)\n")

G = const.G / (1000**3)
v = ((G*m_p)/r_p)**0.5

print(f"Orbital velocity : {v} km/s  \n")


print("3)\n")

p_1 = np.array([-3339.883,-4854.827])
t_1 = 209.808823

p_2 = np.array([-465.010,-5874.345])
t_2 = 209.8056903

t   = 209.8304691

d_1 = (t - t_1) * const.c_km_pr_s
d_2 = (t - t_2) * const.c_km_pr_s

Positions = [p_1,p_2]
Distances = [d_1,d_2]

pos = np.array(TrilaterationAlgorithm(r_p,Positions,Distances))
print(f"x : {pos[0]}, y : {pos[1]}\n")

d_1_t = np.linalg.norm(pos - p_1)
d_2_t = np.linalg.norm(pos - p_2)

#y = (p_1[0]*(h**2 + r_p**2 - d_2**2) - p_2[0]*(h**2 + r_p**2 - d_1**2))/(2*(p_2[1] * p_1[0] - p_1[1] * p_2[0]))
print(f"Dist 1 MCast : {d_1} / Dist 1 Trilateration {d_1_t}")
print(f"Dist 2 MCast : {d_2} / Dist 2 Trilateration {d_2_t}\n")