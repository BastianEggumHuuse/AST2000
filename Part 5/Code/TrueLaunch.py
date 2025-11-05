# Regular imports
import numpy             as np
import matplotlib.pyplot as plt
from numba import njit
import pickle as pkl 

from CalculatingTra import Main as Coast

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils

from ast2000tools.space_mission import SpaceMission

#Unpickling launch
with open("Mission.pkl", 'rb') as file:
    mission = pkl.load(file)


t = mission.time_after_launch
fuelMass = 3535.634180709576

# Simulation variables
R_sim = np.array([mission._position_after_launch])
V_sim = np.array([mission._velocity_after_launch])

# Launch variables
R_launch = np.array([mission._position_after_launch])
V_launch = np.array([mission._velocity_after_launch])

# Simulation parameters
FilePath = "NumericalOrbitData.npz"
npz = np.load(FilePath)
config       = npz["config"]
R_planets = npz["r"]
Info = (float(config[0]), int(config[2]), npz["OrbitTimes"], R_planets)
dt = 1/100000
dT = float(config[1])
M = np.zeros(mission.system._number_of_planets + 1)
M[:-1] = mission.system.masses
M[-1] = mission.system.star_mass

# List of boost times (simulation)
boost_t_sim = np.arange(t,t+0.01, dt)


# List of boost dvs (simulation)
boost_v_sim = np.zeros((len(boost_t_sim),2))
# List of boost times (launch)
boost_t_launch = np.ones(len(boost_t_sim))*dt

# List of boost dvs (launch)
boost_v_launch = np.zeros((len(boost_t_sim),2))

# Beginning travel
travel = mission.begin_interplanetary_travel()
travel.restart()
travel.verbose = False

for i in range(len(boost_t_launch)-1):

    #print("\nNew trajectory step:\n")

    # Orienting ourselves
    t_o,R_o,V_o = travel.orient()

    # if t_o == t:
    #     print('Yes t')
    # else: 
    #     print(f'NEI t: {t_o-t}')
    # if R_o[0] == R_sim[-1][0] and R_o[1] == R_sim[-1][1]:
    #     print('Yes R')
    # else: 
    #     print(f'NEI R: {R_o}, {R_sim[-1]}')
    # if V_o[0] == V_sim[-1][0] and V_o[1] == V_sim[-1][1]:
    #    print('Yes V')
    # else: 
    #     print(f'NEI V: {V_o-mission._velocity_after_launch}, {V_sim[-1]-mission._velocity_after_launch }')

    # Updating launch data
    R_launch = np.concatenate((R_launch,np.array([R_o])))
    V_launch = np.concatenate((V_launch,np.array([V_o])))

    # Boosting
    travel.boost(boost_v_launch[i])
    # Boosting sim
    if(i < len(boost_t_sim) - 1):
        V_sim[-1] += boost_v_sim[i]

    # Coasting
    travel.coast(boost_t_launch[i+1])
    # Coasting sim
    if(i  < len(boost_t_sim) - 1):
        N = round((boost_t_sim[i+1] - boost_t_sim[i])/dt)
        R_c, V_c, A_c, t_c  = Coast(R_sim[-1], V_sim[-1], dt, dT, N, M, t, Info)
        R_sim = np.concatenate((R_sim,R_c))
        V_sim = np.concatenate((V_sim,V_c))
        t += t_c

    

# Orienting ourselves
t_o,R_o,V_o = travel.orient()

# One final update
R_launch = np.concatenate((R_launch,np.array([R_o])))
V_launch = np.concatenate((V_launch,np.array([V_o])))

# Plotting

a_launch_0 = (V_launch[2] - V_launch[1])/dt
print("acc_0", a_launch_0)


fig, ax = plt.subplots()

ax.plot(R_planets[0][0],R_planets[1][0])
ax.plot(R_planets[0][1],R_planets[1][1])

ax.plot(R_sim[:,0],R_sim[:,1],'.', color = 'red')
print(R_sim)
ax.plot(R_launch[:,0],R_launch[:,1], '.', color = 'green')

plt.show()
