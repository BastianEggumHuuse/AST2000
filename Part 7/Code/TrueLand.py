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


Landing.adjust_parachute_area(13)
Landing.look_in_direction_of_planet(1)

t,r,v = Landing.orient()
dv = np.array([0,0,-187]) -594*(np.cross(r,np.array([0,0,1])))/np.linalg.norm(r) - v
print(Landing.orient())
Landing.start_video()
Landing.verbose = True
Landing.launch_lander(dv)
Landing.fall(800)
Landing.deploy_parachute()
Landing.fall(50000)
Landing.finish_video('TrueLand_2.xml')

