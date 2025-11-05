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

#Unpickling launch
with open("Mission.pkl", 'rb') as file:
    mission = pkl.load(file)
FilePath = "NumericalOrbitData.npz"
FIndR = NumericalOrbitFunction(FilePath)
T_0 = mission.time_after_launch
R_1 = FIndR(T_0, 0)
R_2 = FIndR(T_0 + FIndR.dt, 0)
def Lerp(R_0,R_1, I):
    dR_x = -R_0[0] + R_1[0]
    dR_y = -R_0[1] + R_1[1]
    
    return R_0[0] + I* dR_x, R_0[1] + I*dR_y
R_p_x, R_p_y = Lerp(R_1, R_2, (25%10)/(10)) 

plt.plot((R_1[0],R_2[0]),(R_1[1],R_2[1]))
plt.plot(R_p_x,R_p_y, 'o')
plt.show()