# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune

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

def Lerp(R_0,R_1, I):

    # Finding the difference between R_0 and R_1
    dR_x = R_1[0] - R_0[0]
    dR_y = R_1[1] - R_0[1]
    
    # Lerping between R_0 and R_1 for both coordinates
    R_x = R_0[0] + I * dR_x
    R_y = R_0[1] + I * dR_y

    return R_x,R_y 

# Launch variables
R_launch = np.array([mission._position_after_launch])
V_launch = np.array([mission._velocity_after_launch])
t = mission.time_after_launch

# Simulation parameters
FilePath     = "NumericalOrbitData.npz"
FindR        = NumericalOrbitFunction(FilePath)
npz          = np.load(FilePath)
config       = npz["config"]
R_planets    = npz["r"]

M            = np.zeros(mission.system._number_of_planets + 1)
M[:-1]       = mission.system.masses
M[-1]        = mission.system.star_mass


# List of boost times 
boost_t_launch = np.array([
    0,
    0.24,
    0.15,
    0.40,
    0.15,
    0.71
    ])

# List of boost dvs 
boost_v_launch = np.array([
    (0.43,-0.14),
    (0,0),
    (0.1,0),
    (0,0),
    (-0.0855,0),
    (0,0),
    ])

# Instantiating interplanetary travel.
travel = mission.begin_interplanetary_travel()
travel.restart()

for i in range(len(boost_t_launch)-1):
    
    # Orienting ourselves
    t_o,R_o,V_o = travel.orient()
    # Updating launch data
    R_launch = np.concatenate((R_launch,np.array([R_o])))
    V_launch = np.concatenate((V_launch,np.array([V_o])))

    # Boosting
    travel.boost(boost_v_launch[i])
    
    # Coasting
    travel.coast(boost_t_launch[i+1])

# Orienting ourselves a final time
t_o,R_o,V_o = travel.orient()
R_launch = np.concatenate((R_launch,np.array([R_o])))
V_launch = np.concatenate((V_launch,np.array([V_o])))

# Finding position of planet at the end of the coasting
R_0 = FindR(t_o, 1)
R_1 = FindR(t_o+FindR.dt, 1)
I = (t_o % FindR.dt)/FindR.dt 
R_p = Lerp(R_0,R_1, I)


# Initalizing plotting
fig, ax = plt.subplots()
# Plotting the planets trajectories
ax.plot(R_planets[0][0],R_planets[1][0])
ax.plot(R_planets[0][1],R_planets[1][1])
# Plotting the destination planet
ax.plot(R_p[0], R_p[1], 'o', color = 'blue' )
# Plotting the rocket trajectory
ax.plot(R_launch[:,0],R_launch[:,1], '.', color = 'green')

# Plotting the area we need to be in to perform orbital injection maneuver.
l = np.linalg.norm(R_launch[-1][-1])*(M[1]/(M[-1]*10))**0.5
OrbitRadi = plt.Circle(R_p, l, color = 'blue', fill = False, ls = '-')
ax.add_patch(OrbitRadi)
plt.show()


"""
Down here, we perform our orbital injection maneuver.
"""

# Finding distance between planet and rocket
r = R_o - R_p
r_hat = r/np.linalg.norm(r)
# Finding tangential vector
e = -np.array([r_hat[1],-r_hat[0]])

# Computing v_stable
v_stable = ((const.G_sol * mission.system.masses[1]/(np.linalg.norm(r)))**0.5) * e + FindR.GetVelocity(t_o, 1)
dV = (v_stable - V_o) 
# Boosting
travel.boost(dV)

# Orienting one final (final) time.
t_o,R_o,V_o = travel.orient()

# Recording our destination
travel.record_destination(1)

# Re-pickling launch
with open ("Mission.pkl", 'wb') as file:
    pkl.dump(mission, file)

"""
Output:

Note: Existing recorded destination was cleared.
The spacecraft and system were reset to initial state.
Performed automatic orientation:
Time: 2.19991 yr
Position: (-0.209949, -2.7575) AU
Velocity: (5.67544, -3.02612) AU/yr
Spacecraft boosted with delta-v (0.43, -0.14) AU/yr (1802.41 kg of fuel was used).
Spacecraft coasted for 0.24 yr.
Performed automatic orientation:
Time: 2.43991 yr
Position: (1.19449, -2.77097) AU
Velocity: (5.53305, 1.20127) AU/yr
Spacecraft boosted with delta-v (0, 0) AU/yr (0 kg of fuel was used).
Spacecraft coasted for 0.15 yr.
Performed automatic orientation:
Time: 2.58991 yr
Position: (1.9718, -2.49326) AU
Velocity: (4.78958, 2.44729) AU/yr
Spacecraft boosted with delta-v (0.1, 0) AU/yr (290.962 kg of fuel was used).
Spacecraft coasted for 0.4 yr.
Performed automatic orientation:
Time: 2.98991 yr
Position: (3.41537, -1.08291) AU
Velocity: (2.25319, 4.2683) AU/yr
Spacecraft boosted with delta-v (0, 0) AU/yr (0 kg of fuel was used).
Spacecraft coasted for 0.15 yr.
Performed automatic orientation:
Time: 3.13991 yr
Position: (3.67604, -0.423521) AU
Velocity: (1.22635, 4.48826) AU/yr
Spacecraft boosted with delta-v (-0.0855, 0) AU/yr (224.528 kg of fuel was used).
Spacecraft coasted for 0.71 yr.
Performed automatic orientation:
Time: 3.84991 yr
Position: (2.91051, 2.50762) AU
Velocity: (-2.99942, 3.09788) AU/yr
Spacecraft boosted with delta-v (-0.491576, 0.405236) AU/yr (1145.56 kg of fuel was used).
Performed automatic orientation:
Time: 3.84991 yr
Position: (2.91051, 2.50762) AU
Velocity: (-3.49099, 3.50312) AU/yr
Recorded interplanetary travel destination:
Time: 3.84991 yr
Position: (2.91051, 2.50762) AU
Velocity: (-3.49099, 3.50312) AU/yr
"""