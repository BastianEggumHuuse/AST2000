import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import  ast2000tools.constants as const
from ast2000tools.relativity import RelativityExperiments


#1.4.A
M_s = const.m_sun*const.G/(const.c**2)
print(f'Mass sun in m :{M_s}')

Ratio_M2R = M_s/const.R_sun
print(f'ratio M/r for sun : {Ratio_M2R}')

#1.4.C 
M_j = 5.972e24 *const.G/(const.c**2)
print(f'Mass earth in m: {M_j}')

Ratio_Mj2R = M_j/const.R_sun
print(f'ratio M/r for earth : {Ratio_Mj2R}')

#1.4.D

l_j = 500e-9 *(1 + Ratio_M2R) *Ratio_Mj2R
print(f'Endringen i bølgelengden/l_0 som treffer jorda er: {Ratio_Mj2R}')
print(f'Endringen i bølgelengden fra lys på 500nm fra solen som treffer jorden er {l_j*1e9} nm')
