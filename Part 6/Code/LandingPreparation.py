# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune

# Regular imports
import numpy             as np
import matplotlib.pyplot as plt
from numba import njit
import pickle as pkl 

from GeneralizedLaunch import NumericalOrbitFunction

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

def CoordinateAtTime(coordinate_at_zero,elapsed_time,p_theta):

    """
    Function that computes the position of a point on the planet after a given time.

    parameters:
    coordinate_at_zero : array(float) | coordinates (spherical) of the point we wish to find, at time = 0
    elapsed_time       : float        | the time we want to find the point for
    p_theta            : float        | angular velocity of points on the planet (how fast the planet spins)

    returns:
    array(float) | the same point on the planet as seen in coordinate_at_zero, but after the given time has elapsed
    """

    r_0 = coordinate_at_zero[0]
    theta_0 = coordinate_at_zero[1]

    r = r_0
    theta = theta_0 + p_theta * elapsed_time

    coordinate_at_time = np.array((r,theta,0))
    return coordinate_at_time

def AngleAtZero(coordinate_at_time,elapsed_time,p_theta):

    """
    Function that computes the position of a point on the planet at time = 0

    parameters:
    coordinate_at_time : array(float) | coordinates (spherical) of the point we wish to find, at given time
    elapsed_time       : float        | the time which has elapsed
    p_theta            : float        | angular velocity of points on the planet (how fast the planet spins)

    returns:
    float | the angle phi of the point given with coordinate_at_time, but at time = 0
    """

    theta = coordinate_at_time[0]
    
    theta_0 = theta - p_theta * elapsed_time

    return theta_0

def FindAngle(landing_sequence):

    t,r,v = landing_sequence.orient()

    # Finding angle
    if(r[0] == 0):
        # Not included in the flowchart.
        # If x-coordinate is 0, this division is illegal, so we introduce a small number instead
        theta = np.arctan(r[1]/(0.000001))
    else:
        theta = np.arctan(r[1]/r[0])
    if(r[0] < 0):
        theta += (r[1]/abs(r[1]))*np.pi 

    return theta,t

if __name__ == "__main__":

    #Unpickling launch
    with open("Mission.pkl", 'rb') as file:
        mission = pkl.load(file)

    # Beginning landing sequence
    landing = mission.begin_landing_sequence()

    # Getting initial conditions
    t_0,r_vec_0,v_vec_0 = landing.orient()