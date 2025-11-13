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
from GeneralizedLaunch import NumericalOrbitFunction

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
dt = 1/1000000
dT = float(config[1])
M = np.zeros(mission.system._number_of_planets + 1)
M[:-1] = mission.system.masses
M[-1] = mission.system.star_mass
FindR = NumericalOrbitFunction(FilePath)

# List of boost times (simulation)
boost_t_sim = np.array([
    0,
    0.24,
    0.47,
    0.479,
    ])

# List of boost dvs (simulation)
boost_v_sim = np.zeros((len(boost_t_sim),2))
# List of boost times (launch)
boost_t_launch = np.array([
    0,
    0.24,
    0.15,
    0.40,
    0.15,
    0.71
        ])

# List of boost dvs (launch)
boost_v_launch = np.array([
    (0.43,-0.14),
    (0,0),
    (0.1,0),
    (0,0),
    (-0.0855,0),
   
    (0,0),
    (0,0),
    (0,0),
    (0,0),
    (0,0)
     ])
# Beginning travel
travel = mission.begin_interplanetary_travel()
travel.restart()


for i in range(len(boost_t_launch)-1):

    #print("\nNew trajectory step:\n")
    
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
R_p = FindR(t_o, 1)
ax.plot(R_p[0], R_p[1], 'o', color = 'blue' )
ax.plot(R_sim[:,0],R_sim[:,1], color = 'red')

ax.plot(R_launch[:,0],R_launch[:,1], '.', color = 'green')

l = np.linalg.norm(R_launch[-1][-1])*(M[1]/(M[-1]*10))**0.5
OrbitRadi = plt.Circle(R_p, l, color = 'blue', fill = False, ls = '-')
ax.add_patch(OrbitRadi)


r = R_o - R_p

r_hat = r/np.linalg.norm(r)
e = -np.array([r_hat[1],-r_hat[0]])

v_stable = ((const.G_sol * mission.system.masses[1]/(np.linalg.norm(r)))**0.5) * e + FindR.GetVelocity(t_o, 1)
dV = (v_stable - V_o) 
print(dV)
travel.boost(dV)

t_o,R_o,V_o = travel.orient()

travel.record_destination(1)

with open ("Mission.pkl", 'wb') as file:
    pkl.dump(mission, file)