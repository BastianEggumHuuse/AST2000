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

#Doing everything else 
m_p = mission.system.masses[1] * const.m_sun
G   = const.G
#intializing lanidng
landing = mission.begin_landing_sequence()
landing.verbose = False
#orienting for intial possistion
t,r_0,v_0 = landing.orient()




#starting video
landing.look_in_direction_of_motion()
landing.start_video()
#setting parameters 
N = 1000
R = []
r_ø = 2.7e6
print(f'Ønsket distangse over bakken : {r_ø/1000 - landing.system.radii[1]}km' )
#a loop where we slowly fall down bossting with a small bost then faling some more till we are close enugh
for n in range(N):


    #orienting and falling
    landing.fall(60)
    t,r,v = landing.orient()
    #cheking if we are as close as we want
    if np.linalg.norm(r) < r_ø:
        v_s = np.sqrt((const.G*mission.system.masses[1]*const.m_sun/np.linalg.norm((r))))
        v_hat = -np.array((r[1], -r[0],0))/np.linalg.norm(r)
        dv = v_s*v_hat - v
        
        landing.boost(dv)
       
        break
    
    #Apending to list
    R.append(r)
    v_hat = v / np.linalg.norm(v_0)
    #boosting 
    dv = -1000 * v_hat
    dv_i = dv/N
    landing.boost(dv_i)
#I dont like lists ok
R = np.array(R)
#falling to get the video better, also this then becomes the basis for the rest
landing.fall(30000)
landing.finish_video()


#pickling here 
with open ("Landing.pkl", 'wb') as file:
    pkl.dump(landing, file)
#creating test arrays
R_T = np.zeros(100)
T_T = np.zeros(100)
#Looping over and storing mean distances for for a timperiod for plotting
for i in range (100):
    r_mean = 0
    for k in range(100):
       
        landing.fall(20)
        
        t_o,r,v = landing.orient()
        r_mean += np.linalg.norm(r)
    #taking the mean, by dividing by the number of steps
    R_T[i] = r_mean/100
    
    T_T[i] = t_o


#plotng
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
"""
kjøretidseksempel
Note: Ongoing interplanetary travel was terminated.
Ønsket distangse over bakken : 395.4056984029403km
XML file landing_video.xml was saved in XMLs/.
It can be viewed in MCAst.
"""