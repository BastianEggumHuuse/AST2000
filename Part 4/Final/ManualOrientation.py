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

"""
Code that runs this program: python ManualOrientation.py

Output:

Initializing motor...
Calculated Force per motor : 1.26686e-10, Calculated Fuel Consumption per motor : 2.95044e-14
Calculated Force           : 2.11143e+06, Calculated Fuel Consumption      : 4.91740e+02

Generalized Position in solar system frame at t = 1: [x : -1.302 AU, y : 2.42e+00 AU]
Generalized Velocity in solar system frame at t = 1: [x : -6.414 AU/Y, y : -0.728 AU/Y]

Specialized Position in solar system frame at t = 1: [x : 2.814 AU, y : 6.98e-05 AU]
Specialized Velocity in solar system frame at t = 1: [x : 2.480 AU/Y, y : 5.856 AU/Y]

Rocket was moved down by 6524.35 m to stand on planet surface.
New launch parameters set.
Launch completed, reached escape velocity in 379.294 s.
Your spacecraft position was satisfyingly calculated. Well done!
*** Achievement unlocked: No free launch! ***
Picture written to Sky_Image.png.
Pointing angle after launch correctly calculated. Well done!
Velocity after launch correctly calculated. Well done!
Position after launch correctly calculated. Well done!
Your manually inferred orientation was satisfyingly calculated. Well done!
*** Achievement unlocked: Well-oriented! ***
"""