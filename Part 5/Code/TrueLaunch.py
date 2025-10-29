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
boost_t_sim = np.array([
    t,
    t + 0.24,
    t + 0.47,
    t + 0.479
    ])

# List of boost dvs (simulation)
boost_v_sim = np.array([
    (0,0),
    (0.14,0.29),
    (-0.1,0),
    (0,0)
    ])

# List of boost times (launch)
boost_t_launch = np.array([
    t,
    t + 0.0001,
    t + 0.24,
    t + 0.47,
    #t + 0.479,
    t + 1,
    ])

# List of boost dvs (launch)
boost_v_launch = np.array([
    #(0.4,-0.7),
    (0.0,0.0),
    (0.0,0.0),
    (0.14,0.29),
    (-0.1,0),
    (0,0),
    ])

# Beginning travel
travel = mission.begin_interplanetary_travel()

for i in range(len(boost_t_launch)-1):

    print("\nNew trajectory step:\n")

    # Orienting ourselves
    t_o,R_o,V_o = travel.orient()

    # Updating launch data
    R_launch = np.concatenate((R_launch,np.array([R_o])))
    V_launch = np.concatenate((V_launch,np.array([V_o])))

    # Boosting
    travel.boost(boost_v_launch[i])
    # Boosting sim
    if(i < len(boost_t_sim) - 1):
        V_sim[-1] += boost_v_sim[i]

    # Coasting
    travel.coast_until_time(boost_t_launch[i+1])
    # Coasting sim
    if(i  < len(boost_t_sim) - 1):
        N = int(np.floor((boost_t_sim[i+1] - boost_t_sim[i])/dt))
        R_c, V_c, A_c, t_c  = Coast(R_sim[-1], V_sim[-1], dt, dT, N, M, t, Info)
        R_sim = np.concatenate((R_sim,R_c))
        V_sim = np.concatenate((V_sim,V_c))
        t += t_c

    print()

# Orienting ourselves
t_o,R_o,V_o = travel.orient()

# One final update
R_launch = np.concatenate((R_launch,np.array([R_o])))
V_launch = np.concatenate((V_launch,np.array([V_o])))

# Plotting

fig, ax = plt.subplots()

ax.plot(R_planets[0][0],R_planets[1][0])
ax.plot(R_planets[0][1],R_planets[1][1])

ax.plot(R_sim[:,0],R_sim[:,1],ls = "-")
ax.plot(R_launch[:,0],R_launch[:,1])

plt.show()
