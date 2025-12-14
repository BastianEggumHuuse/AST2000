import numpy as np
import matplotlib.pyplot as plt

import ast2000tools.utils as utils
from   ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as const
from   ast2000tools.relativity import RelativityExperiments

LorentzFactor = lambda v : 1/(1-v**2)

M = 1
v_shell = 0.993
theta = np.deg2rad(167)

R = 20*M
L_m = R * LorentzFactor(v_shell) * v_shell * np.sin(theta)

def Veff(r):

    return ((1-(2*M)/r) * (1 + (L_m**2)/(r**2)))**0.5

r_E_min = ((L_m**2)/(2*M))*(1+(1-12*(M**2)/(L_m**2))**0.5)
r_E_max = ((L_m**2)/(2*M))*(1-(1-12*(M**2)/(L_m**2))**0.5)

r = np.linspace(2*M,300*M,10000)

r_Rocket   = 20 * M
E_m_Rocket = ((1 - (2*M)/r_Rocket)**0.5) * LorentzFactor(v_shell)

plt.axvline(r_E_max,linestyle = "--",color = "green",label = r"$r_{crit}$")
#plt.axhline(Veff(r_E_min),linestyle = "--",color = "red")
plt.axhline(E_m_Rocket,linestyle = "dotted",color = "royalblue",label = r"$\frac{E}{m}$")

plt.plot(r,Veff(r),color = "darkorange")

plt.xlabel("Distanse fra sentrum \nav det sorte hullet [M]")
plt.ylabel(r"$V_{eff}$")

plt.tight_layout()
plt.legend()
plt.show()