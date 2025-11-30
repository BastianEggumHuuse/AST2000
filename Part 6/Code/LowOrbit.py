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
R = []
r_ø = 2.7e6
print(f'Ønsket distangse over bakken : {r_ø/1000 - landing.system.radii[1]}km' )
for n in range(N):


    
    landing.fall(60)
    t,r,v = landing.orient()
    
    if np.linalg.norm(r) < r_ø:
        v_s = np.sqrt((const.G*mission.system.masses[1]*const.m_sun/np.linalg.norm((r))))
        v_hat = -np.array((r[1], -r[0],0))/np.linalg.norm(r)
        dv = v_s*v_hat - v
        
        landing.boost(dv)
       
        break
    
    
    R.append(r)
    v_hat = v / np.linalg.norm(v_0)
    dv = -1000 * v_hat
    dv_i = dv/N
    landing.boost(dv_i)
R = np.array(R)
landing.fall(30000)
landing.finish_video()



with open ("Landing.pkl", 'wb') as file:
    pkl.dump(landing, file)
R_T = np.zeros(100)
T_T = np.zeros(100)

for i in range (100):
    r_mean = 0
    for k in range(100):
       
        landing.fall(20)
        
        t_o,r,v = landing.orient()
        r_mean += np.linalg.norm(r)
    
    R_T[i] = r_mean/100
    
    T_T[i] = t_o



plt.plot(0,0,'o')
plt.xlabel('X [m]')
plt.ylabel('Y [m]')
plt.plot(R[:,0], R[:,1])
plt.show()
plt.plot(T_T,R_T)
plt.ylabel('Gjennomsnitt radius [m]')
plt.xlabel('TId [s]')
plt.show()
landing.orient()
print(np.linalg.norm(r))