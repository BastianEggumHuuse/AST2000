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

with open("Landing.pkl", 'rb') as file:
    Landing = pkl.load(file)

Landing.fall(1230)

Landing.adjust_parachute_area(13)
Landing.look_in_direction_of_motion()

t,r,v = Landing.orient()
dv = v_0 = np.array([0,0,-180]) -870*(np.cross(r,np.array([0,0,1])))/np.linalg.norm(r) - v
print(dv)
Landing.look_in_direction_of_motion()
print(Landing.orient())
Landing.start_video()
Landing.verbose = True
Landing.launch_lander(dv)
Landing.fall(800)
Landing.look_in_direction_of_motion()
Landing.deploy_parachute()
Landing.fall(5000)
Landing.look_in_direction_of_motion()
Landing.finish_video('TrueLand_1.xml', 10000)

