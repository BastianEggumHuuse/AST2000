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
r_h = np.linalg.norm(r_vec_1)
print(f"Height : {h} km\n")

print("2)\n")

G = const.G / (1000**3)
v = ((G*m_p)/r_p)**0.5

print(f"Orbital velocity : {v} km/s  \n")


print("3)\n")

p_1 = np.array([-3339.883,-4854.827])
t_1 = 209.8088823

p_2 = np.array([-465.010,-5874.345])
t_2 = 209.8056903

t   = 209.8304691

def findlocation(p_1,p_2,t_0,t_1,t_2, kvandrant3 = True):
    d_1 = (t_0 - t_1) * const.c_km_pr_s
    d_2 = (t_0 - t_2) * const.c_km_pr_s

    
    phi_s1 = np.arctan(p_1[1]/p_1[0]) 
    phi_p1 = np.arccos(-(d_1**2-np.linalg.norm(p_1)**2-r_p**2)/(2*np.linalg.norm(p_1)*r_p))

    phi_s2 = np.arctan(p_2[1]/p_2[0]) 
    phi_p2 = np.arccos(-(d_2**2-np.linalg.norm(p_2)**2-r_p**2)/(2*np.linalg.norm(p_2)*r_p))
    
    if kvandrant3:
        phi_s1 += np.pi
        phi_s2 += np.pi


    theta_11 = phi_s1 + phi_p1
    theta_12 = phi_s1 - phi_p1

    theta_21 = phi_s2 + phi_p2
    theta_22 = phi_s2 - phi_p2
    
    print(f'Ship 1 : {theta_11 :3f} og {theta_12 :3f}')
    print(f'Ship 2 : {theta_21 :3f} og {theta_22 :3f}\n')

    dif = 1e-3
    if abs(theta_11 - theta_21) < dif or abs(theta_11 - theta_22) < dif:
        theta = theta_11
    elif abs(theta_12 - theta_21) < dif or abs(theta_12 - theta_22) < dif:
        theta = theta_12
    else:
        raise ValueError('Ingen theta like')


    x = np.cos(theta)*r_p ; y = np.sin(theta)*r_p
    print(f"x : {x}, y : {y}\n")
    pos = np.array([x,y])
    return pos, d_1, d_2
pos, d_1, d_2 = findlocation(p_1, p_2, t, t_1,t_2)
d_1_t = np.linalg.norm(pos - p_1)
d_2_t = np.linalg.norm(pos - p_2)


print(f"Dist 1 MCast : {d_1} / Dist 1 Trilateration {d_1_t}")
print(f"Dist 2 MCast : {d_2} / Dist 2 Trilateration {d_2_t}\n")

print("4)\n")
#converting everything to rel units
def GetRelT(t_1, t_2):
    v_rel = v/const.c_km_pr_s
    t_1_m = t_1 #* const.c_km_pr_s 
    t_2_m = t_2 #* const.c_km_pr_s 
    m_p_r = m_p * G/(const.c_km_pr_s**2) 
    print(m_p_r)

    t_r_1 = t_1_m * ((1-2*m_p_r/(r_p) )/(1-2*(m_p_r)/r_h - v_rel**2))**0.5
    t_r_2 = t_2_m * ((1-2*m_p_r/(r_p) )/(1-2*(m_p_r)/r_h - v_rel**2))**0.5
    
    t_rs_1 = t_r_1#/const.c_km_pr_s
    t_rs_2 = t_r_2#/const.c_km_pr_s
    return t_rs_1, t_rs_2
t_r1, t_r2 = GetRelT(t_1,t_2)
print(f'Tiden i sekunder i bakke refferanse system sattelitt 1 = {t_r1}')
print(f'Tiden i sekunder i bakke refferanse system sattelitt 2 = {t_r2}\n')


print("5)\n")

pos_r, d_1_r, d_2_r = findlocation(p_1, p_2, t, t_r1,t_r2)

print(f'Pos Med Relativstisk Teori:   {pos_r}')
print(f'Pos uten Relativistisk Teori: {pos}')
print(f'Differanse =                        {np.linalg.norm(pos_r - pos)*1000}m ')

print("6)\n")


p_1 = np.array([4200.968,4132.315])
t_1 = 13075.0391380
p_2 = np.array( [5704.303,1478.206])
t_2 = 13075.0360543
t = 13075.0611044

pos, d_1, d_2 = findlocation(p_1, p_2, t, t_1,t_2, kvandrant3=False)

t_r1, t_r2 = GetRelT(t_1,t_2)

print(f'Tiden i sekunder i bakke refferanse system sattelitt 1 = {t_r1}')
print(f'Tiden i sekunder i bakke refferanse system sattelitt 2 = {t_r2}\n')


pos_r_1, d_1_r, d_2_r = findlocation(p_1, p_2, t, t_r1,t_r2,kvandrant3=False)
print(f'Pos Med Relativstisk Teori:   {pos_r_1}')
print(f'Pos uten Relativistisk Teori: {pos}')
print(f'Differanse =                        {np.linalg.norm(pos_r_1 - pos)*1000}m ')
print(f'før - etter:  {np.linalg.norm(pos_r-pos_r_1)*1000}')