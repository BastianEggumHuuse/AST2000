import numpy as np
import matplotlib.pyplot as plt

import ast2000tools.utils as utils
from   ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as const
import scipy.constants as const2
from   ast2000tools.relativity import RelativityExperiments

m   = 1e6
v_a = 0.177024
N   = 5.39693e41

c = const.c
h = const2.h

f = (2*m*c**2)/((1-v_a**2)**0.5 * N * h)

l = c / f

print(f"lambda         : {l * 1e9:.5f} nm")

dl = l*(np.sqrt((1+v_a)/(1-v_a))-1)

print(f"delta lambda   : {dl * 1e9:.5f} nm")
print(f"shifted lambda : {(l - dl) * 1e9:.5f} nm")