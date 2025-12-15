import numpy as np
import matplotlib.pyplot as plt

import ast2000tools.utils as utils
from   ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as const
from   ast2000tools.relativity import RelativityExperiments
LorentzFactor   = lambda v : 1/(1-v**2)
v_shell = 0.993
theta = np.deg2rad(167)
M = 1
L_m = 20*M*LorentzFactor(v_shell)*v_shell * np.sin(theta)

def Veff(r):
    return ((1 - (2*M)/r)*(1 + (L_m)**2/(r**2)))**0.5
print(L_m)
Min_r =  1/2 * L_m**2/M *(1 + (1 - 12*M**2/L_m**2)**0.5)
Maks_r = 1/2 * L_m**2/M *(1 - (1 - 12*M**2/L_m**2)**0.5)
print(Veff(Min_r))

r = np.linspace(2, 2*102544,1000000)*M
E = Veff(r)
print(Min_r)
plt.plot(r,E)
#plt.axvline(Min_r)
#plt.axhline(Veff(Min_r))
plt.xlabel('R [M]')
plt.ylabel('Energi []')
plt.show()