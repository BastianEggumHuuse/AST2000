# Regular imports
import numpy             as np
import matplotlib.pyplot as plt
from numba import njit
import pickle as pkl 

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils

from ast2000tools.space_mission import SpaceMission
from GeneralizedLaunch import NumericalOrbitFunction

#Unpickling launch
with open("Mission.pkl", 'rb') as file:
    mission = pkl.load(file)

m_p = mission.system.masses[1] * const.m_sun
G   = const.G

landing = mission.begin_landing_sequence()
landing.verbose = False

t,r_0,v_0 = landing.orient()
v_hat = v_0 / np.linalg.norm(v_0)
dv = -200 * v_hat

N = 10
dv_i = dv/N



N = 10
r_mean = []
for n in range(N):
    print(f"Running generation {n}")

    I = 10
    r_mean_i = 0
    for i in range(I):
        landing.fall(3600*24 * 3.5)
        t,r,v = landing.orient()
        r_mean_i += np.linalg.norm(r)

    r_mean_i *= 1/I

    r_mean.append(r_mean_i)
    print(f"Mean for gen {n} is {r_mean_i}")
    landing.boost(dv_i)
    
for r in r_mean:

    print(f"mean r : {r:10.3f}")


plt.plot(np.arange(len(r_mean)),np.ones(len(r_mean))*np.linalg.norm(r_0))
plt.plot(np.arange(len(r_mean)),r_mean)
plt.show()