# Ikke brukt kodemal!!
# Skrevet av Bastian Eggum Huuse og Bendik Thune

import sys
import numpy as np
from numba import njit
from PIL import Image

import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

from GeneralizedLaunch import main
from ImageAnalysis import CompareImageRange
from DopplerShift import Main as DopplerVelocity
from Trilateration import TrilaterationAlgorithm

if __name__ ==  "__main__":

    # Ast init
    seed = utils.get_seed('bmthune')
    mission = SpaceMission(seed)

    # Defining launch time
    t_0 = 1

    # Running Launch
    t_1 = main(mission,t_0)


    # Collecting information:
    # Rotation:
    InputPath = "Sky_Image.png"
    mission.take_picture(InputPath,"himmelkule.npy")
    InputImage = Image.open(InputPath)
    InputImageArray = np.array(InputImage)
    SkyImageData = np.load("SkyImageData.npy")
    # Velocity:
    lambda_0 = 656.3 #Hydrogen spectral line
    lambda_1, lambda_2 = mission.star_doppler_shifts_at_sun
    lambda_3, lambda_4 = mission.measure_star_doppler_shifts()
    dlambda = (lambda_1,lambda_2,lambda_3,lambda_4)
    phi_1, phi_2 = mission.star_direction_angles
    phi_1, phi_2 = (np.deg2rad(phi_1), np.deg2rad(phi_2))
    # Position
    Distances = mission.measure_distances()

    # Performing computations
    Orientation_Rotation = CompareImageRange(InputImageArray,SkyImageData,N = 359)
    Orientation_Velocity = DopplerVelocity(dlambda, lambda_0, phi_1, phi_2)
    Orientation_Position = TrilaterationAlgorithm(t_1,Distances)

    # Verifying results
    mission.verify_manual_orientation(Orientation_Position,Orientation_Velocity,Orientation_Rotation)