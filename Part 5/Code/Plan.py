# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune

# Regular imports
import numpy             as np
import matplotlib.pyplot as plt
from numba import njit
import pickle as pkl 

from GeneralizedLaunch import NumericalOrbitFunction,FuelRocket
from GeneralizedLaunch import NumericalOrbitFunction
from CalculatingTra import Main as Coast

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

# Importing the numerical planet positions and velocities
FilePath = "NumericalOrbitData.npz"
PlanetPositionFunction = NumericalOrbitFunction(FilePath)

# We want to use numba to speed up the computing process, but 
# NumericalOrbitFunction is a class, which numba isn't a huge fan of
# Therefore, we also load the data locally, to find planet positions.
npz = np.load(FilePath)
TotalTime    = float(npz["config"][0])
NumSteps     = int(npz["config"][2])
OrbitTimes   = npz["OrbitTimes"]
r            = npz["r"]
# Storing all this data in a handy tuple
Info = (TotalTime, NumSteps, OrbitTimes, r)

"""
This program is not very interesting, so it has no corresponding flowchart.
It's structure is essentially just

1) Do a boost                            (which is one line of code)
2) Calculate fuel lost                   (which we did in part 1)
3) Calculate trajectory until next boost (which we did earlier this part)
4) Repeat
"""

#Unpickling launch
with open("Mission.pkl", 'rb') as file:
    mission = pkl.load(file)

def FuelEstimate(fuel_mass, mission, dV):
    """
    Function that calculates how much fuel is neccesary to boost the rocket.

    Parameters:
    fuel_mass : float        | Mass of remaining rocket fuel
    mission   : SpaceMission | Our SpaceMission instance
    dV        : Array(float) | The boost we wish to perform

    returns:
    float | How much fuel (in kg) is needed to perform the boost dV
    """ 

    # These are the same parameters our rocket uses
    NumMotors = int((1000000**3)/60) # 1/10 qube meter grid :)
    NumParticles = 10**5

    # Finding the length of the boost, and converting to m/s
    dv = np.linalg.norm(dV) * (const.AU / (60*60*24*365))

    # Running a fuelrocket instance (written in part 1)
    f = FuelRocket(fuel_mass,dv,NumMotors,mission,NumParticles)
    f.TimeLoop()

    # Returning fuel diff.
    return fuel_mass - f.FuelMass

"""
This part of the code simply visualises the different positions of our home planet and destination planet
over time. We do this to find a satisfactory launch time.
"""

t = np.arange(0,5, 0.2)

# Importing the numerical planet positions and velocities
FilePath = "NumericalOrbitData.npz"
Numericalsim = NumericalOrbitFunction(FilePath)
R_p = Numericalsim.range(0,10)

# Plotting the movement of both planets
plt.plot(R_p[0][0],R_p[1][0])
plt.plot(R_p[0][1],R_p[1][1])

# Looping over all times and plotting position of home and destination planets,
# along with a line between them, connecting them.
colors = ['red','orange', 'yellow', 'green', 'blue', 'indigo', 'purple', 'brown']
for i in range(len(t)):
    
    m = i//len(colors)

    color = colors[i-m*len(colors)]

    R_0 = Numericalsim(t[i], 0) #planet 0, planet 1
    R_1 = Numericalsim(t[i], 1)
    plt.plot((R_0[0],R_1[0]),(R_0[1],R_1[1]), color =color)
    plt.plot(R_0[0], R_0[1], 'o',color =color)
    plt.plot(R_1[0], R_1[1], 'o',color =color)
    plt.xlabel('Distanse langs x-aksen [AU]')
    plt.ylabel('Distanse langs y-aksen [AU]')
 

"""
From here on out comes the actual program. :)
"""

# Launch time we landed on
t = 11*0.2

# Getting initial positions for poth planets
R_0 = Numericalsim(t, 0)
R_1 = Numericalsim(t, 1)

# Setting simulation parameters
dt = 1/1000000
dT = Numericalsim.dt
M = np.zeros(mission.system._number_of_planets + 1)
M[:-1] = mission.system.masses
M[-1] = mission.system.star_mass

# Initializing arrays
R = np.array([mission._position_after_launch])
V = np.array([mission._velocity_after_launch])
t = mission.time_after_launch

# Getting position of planets after launch
R_0_s = Numericalsim(t, 0)
R_1_s = Numericalsim(t, 1)

# Setting initial fuel mass
fuelMass = 5000

# Defining set of boosts, along with the time t they occur at.
DV = [
    (t, np.array([0.4,0])),
    (t + 0.15, np.array([0.3,0.4])),
    (t + 0.39, np.array([0,0.6])),
    (t + 0.79, np.array([0,0])),
    (t + 0.79, np.array([0,0])),
    (t + 0.883, np.array([0,0])),
    ]   

# Looping over all boosts and applying them
for i in range(len(DV)-1):
    # Boosting
    V[-1] += DV[i][1]
    # Coasting
    N = int(np.floor((DV[i+1][0] - DV[i][0])/dt))
    R_c, V_c, A_c, t_c  = Coast(R[-1], V[-1], dt, dT, N, M, t, Info)
    R = np.concatenate((R,R_c))
    V = np.concatenate((V,V_c))
    # Updating time
    t += t_c
    # Updating fuel
    fuelMass -= FuelEstimate(fuelMass, mission, DV[i][1])

# Finding position of planets after simulation
print(f'Remaning fuel after launch :{fuelMass}')

R_0 = Numericalsim(t, 0)
R_1 = Numericalsim(t, 1)

fig, ax = plt.subplots()
ax.plot(R_p[0][0],R_p[1][0])
ax.plot(R_p[0][1],R_p[1][1])
ax.plot(R_0[0], R_0[1], 'o')
ax.plot(R_1[0], R_1[1], 'o')
ax.plot(R_0_s[0], R_0_s[1], 'o')
ax.plot(R_1_s[0], R_1_s[1], 'o')
ax.set_xlabel('x :[AU]')
ax.set_ylabel('y :[AU]')

# Computing distance we have to be from destination planet to perform orbital injection maneuver
l = np.linalg.norm(R[-1])*(M[1]/(M[-1]*10))**0.5
OrbitRadi = plt.Circle(R_1, l, color = 'Lime', fill = False, ls = '-')
ax.add_patch(OrbitRadi)

ax.plot(R_0[0], R_0[1], 'o')
ax.plot(R_1[0], R_1[1], 'o')
ax.plot(R[:,0],R[:,1])   
print(f'Distance to target planet over target distance = {np.linalg.norm(R[-1]-R_1)/l}')
plt.show()

"""
Output:

Remaning fuel after launch :62.43865999963532
Distance to target planet over target distance = 8.481437670001904
"""