import numpy as np
import matplotlib.pyplot as plt

import ast2000tools.utils as utils
from   ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as const
from   ast2000tools.relativity import RelativityExperiments

# Data gathered from MCast
m_p = 1.67262158e-27    
m_n = 1.67492747e-27
m_e = 9.10938188e-31

print(f"Mass of proton {m_p:.3e} kg")
print(f"Mass of neutron {m_n:.3e} kg")
print(f"Mass of electron {m_e:.3e} kg\n")

g_p = (m_n**2 + m_p**2 - m_e**2)/(2*m_p*m_n)
g_e = (m_n**2 + m_e**2 - m_p**2)/(2*m_e*m_n)

v_p = -np.sqrt(1-((2*m_p*m_n)/(m_n**2 + m_p**2 - m_e**2))**2)
v_e = np.sqrt(1-((2*m_e*m_n)/(m_n**2 + m_e**2 - m_p**2))**2)

print(f"Speed of proton after neutron decay {v_p}")
print(f"Speed of electron after neutron decay {v_e}\n")

E_p = g_p * m_p
E_e = g_e * m_e

P_p = g_p * m_p * v_p
P_e = g_e * m_e * v_e

# Now we change systems :)

v_rel = 0.8930
g_rel = 1/np.sqrt(1-v_rel**2)

E_p_l = g_rel * E_p + v_rel * g_rel * P_p
E_e_l = g_rel * E_e + v_rel * g_rel * P_e
P_p_l = g_rel * P_p + v_rel * g_rel * E_p
P_e_l = g_rel * P_e + v_rel * g_rel * E_e

print(f"E_p : {E_p_l} kg")
print(f"P_p : {P_p_l} kg\n")

print(f"E_e : {E_e_l} kg")
print(f"P_e : {P_e_l} kg\n")

v_p_l = P_p_l/E_p_l
v_e_l = P_e_l/E_e_l

print(f"v_p_lab         : {v_p_l}")
print(f"v_e_lab         : {v_e_l}")
print(f"Size difference : {(v_e_l - v_rel) / abs(v_p_l-v_rel)}\n")

v_p_l_l = (v_p + v_rel)/(1+v_rel * v_p)
v_e_l_l = (v_e + v_rel)/(1+v_rel * v_e)

print(f"v_p_lab         : {v_p_l_l}")
print(f"v_e_lab         : {v_e_l_l}")
print(f"Size difference : {(v_e_l_l - v_rel) / abs(v_p_l_l-v_rel)}\n")