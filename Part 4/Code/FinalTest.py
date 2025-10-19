# Ikke brukt kodemal!!
# Skrevet av Bastian Eggum Huuse og Bendik Thune

import sys
import numpy as np
from numba import njit

import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

from GeneralizedLaunch import main
from ImageAnalysis import CompareImageRange
from Trilateration import TrilaterationAlgorithm

if __name__ ==  "__main__":

    # Ast init
    seed = utils.get_seed('bmthune')
    mission = SpaceMission(seed)

    # Defining launch time
    t_0 = 1

    # Running Launch
    t_1 = main(mission,t_0)


    # Collecting information
    Image = mission.take_image("Sky_Image.png","himmelkule.npy")
    # Velocity bit goes here
    Distances = mission.measure_distances()

    # Performing computations
    Orientation_Rotation = CompareImageRange(Image,"himmelkule.npy",N = 359)
    # Velocity bit goes here
    Orientation_Position = TrilaterationAlgorithm(t_1,Distances)

    

    
    