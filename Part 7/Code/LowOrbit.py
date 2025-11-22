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





landing.look_in_direction_of_motion()
landing.start_video()

N = 1000
r_mean = []
r_ø = 2.7e6
for n in range(N):

    I = 1
    r_mean_i = 0
    for i in range(I):
        landing.fall(60)
        t,r,v = landing.orient()
        r_mean_i += np.linalg.norm(r)
    if np.linalg.norm(r) < r_ø:
        v_s = np.sqrt((const.G*mission.system.masses[1]*const.m_sun/np.linalg.norm((r))))
        v_hat = -np.array((r[1], -r[0],0))/np.linalg.norm(r)
        dv = v_s*v_hat - v
        print(dv)
        landing.boost(dv)
       
        break
    r_mean_i *= 1/I

    r_mean.append(r_mean_i)
    v_hat = v / np.linalg.norm(v_0)
    dv = -1000 * v_hat
    dv_i = dv/N
    landing.boost(dv_i)

landing.fall(30000)
landing.finish_video()



with open ("Landing.pkl", 'wb') as file:
    pkl.dump(landing, file)
print(landing.orient())
plt.plot(np.arange(len(r_mean)),np.ones(len(r_mean))*np.linalg.norm(r_0))
plt.plot(np.arange(len(r_mean)),r_mean)
plt.show()