# Regular imports
import numpy             as np
import matplotlib.pyplot as plt
from numba import njit
import pickle as pkl 

from GeneralizedLaunch import NumericalOrbitFunction
from CalculatingTra import Main as Coast
# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils

from ast2000tools.space_mission import SpaceMission
FilePath = "NumericalOrbitData.npz"
npz = np.load(FilePath)
config       = npz["config"]
TotalTime    = float(config[0])
NumSteps     = int(config[2])
OrbitTimes   = npz["OrbitTimes"]
r = npz["r"]
Info = (TotalTime, NumSteps, OrbitTimes, r)

#Unpickling launch
with open("Mission.pkl", 'rb') as file:
    mission = pkl.load(file)

def FuelEstimate(fuel_mass, mission, dV):

    NumMotors = int((1000000**3)/60) # 1/10 qube meter grid :)
    NumParticles = 10**5

    dv = np.linalg.norm(dV) * (const.AU / (60*60*24*365))

    f = FuelRocket(fuel_mass,dv,NumMotors,mission,NumParticles)
    f.TimeLoop()

    return fuel_mass - f.FuelMass


#First find a good starting time:
t = np.arange(0,5, 0.2)


Numericalsim = NumericalOrbitFunction(FilePath)
R_p = Numericalsim.range(0,10)

plt.plot(R_p[0][0],R_p[1][0])
plt.plot(R_p[0][1],R_p[1][1])


for t in t:
    R_0 = Numericalsim(t, 0)
    R_1 = Numericalsim(t, 1)
    plt.plot((R_0[0],R_1[0]),(R_0[1],R_1[1]))
    plt.plot(R_0[0], R_0[1], 'o')
    plt.plot(R_1[0], R_1[1], 'o')
 

t = 2.8

R_0 = Numericalsim(t, 0)
R_1 = Numericalsim(t, 1)


dt = 1/1000000
dT = Numericalsim.dt
M = np.zeros(mission.system._number_of_planets + 1)
M[:-1] = mission.system.masses
M[-1] = mission.system.star_mass

R = np.array([mission._position_after_launch])
V = np.array([mission._velocity_after_launch])
t = mission.time_after_launch
R_0_s = Numericalsim(t, 0)
R_1_s = Numericalsim(t, 1)

DV = [
    (t, np.array([0,0])),
    (t + 0.24, np.array([0.14,0.29])),
    (t + 0.47, np.array([-0.1,0])),
    (t + 0.479, np.array([0,0]))
      ]

for i in range(len(DV)-1):
    V[-1] += DV[i][1]
    N = int(np.floor((DV[i+1][0] - DV[i][0])/dt))
    R_c, V_c, A_c, t_c  = Coast(R[-1], V[-1], dt, dT, N, M, t, Info)
    R = np.concatenate((R,R_c))
    V = np.concatenate((V,V_c))
    t += t_c
R_0 = Numericalsim(t, 0)
R_1 = Numericalsim(t, 1)
l = np.linalg.norm(R[-1])*(M[1]/(M[-1]*10))**0.5
fig, ax = plt.subplots()



ax.plot(R_p[0][0],R_p[1][0])
ax.plot(R_p[0][1],R_p[1][1])
ax.plot(R_0[0], R_0[1], 'o')
ax.plot(R_1[0], R_1[1], 'o')
ax.plot(R_0_s[0], R_0_s[1], 'o')
ax.plot(R_1_s[0], R_1_s[1], 'o')

OrbitRadi = plt.Circle(R_1, l, color = 'Lime', fill = False, ls = '-')
ax.add_patch(OrbitRadi)

ax.plot(R_0[0], R_0[1], 'o')
ax.plot(R_1[0], R_1[1], 'o')
ax.plot(R[:,0],R[:,1])   

plt.show()

