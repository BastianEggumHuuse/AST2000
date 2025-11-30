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

def CoordinateAtTime(coordinate_at_zero,elapsed_time,p_phi):

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
    phi_0 = coordinate_at_zero[1]
    theta_0 = coordinate_at_zero[2]

    r = r_0
    phi = phi_0 + p_phi * elapsed_time
    theta = theta_0

    coordinate_at_time = np.array((r,phi,theta))
    return coordinate_at_time

def AngleAtZero(coordinate_at_time,elapsed_time,p_phi):

    """
    Function that computes the position of a point on the planet at time = 0

    parameters:
    coordinate_at_time : array(float) | coordinates (spherical) of the point we wish to find, at given time
    elapsed_time       : float        | the time which has elapsed
    p_theta            : float        | angular velocity of points on the planet (how fast the planet spins)

    returns:
    float | the angle phi of the point given with coordinate_at_time, but at time = 0
    """

    phi = coordinate_at_time[1]
    
    phi_0 = phi - p_phi * elapsed_time

    return phi_0

def FindAngle(landing_sequence):

    t,r,v = landing_sequence.orient()

    # Finding angle
    if(r[0] == 0):
        # Not included in the flowchart.
        # If x-coordinate is 0, this division is illegal, so we introduce a small number instead
        phi = np.arctan(r[1]/(0.000001))
    else:
        phi = np.arctan(r[1]/r[0])
    if(r[0] < 0):
        phi += (r[1]/abs(r[1]))*np.pi 

    while(phi < 0):
        phi += np.pi * 2

    return r,phi,t

if __name__ == "__main__":

    #Unpickling launch
    with open("Mission.pkl", 'rb') as file:
        mission = pkl.load(file)
    
    with open("Landing.pkl", 'rb') as file:
        landing = pkl.load(file)

    # Beginning landing sequence
    landing.look_in_direction_of_planet(1)

    # Getting planet spin
    p_phi = (2*np.pi)/(mission.system.rotational_periods[1]*24*60*60)
    print("Planet Rotational velocity : ",p_phi)

    # Finding the time, position, and velocity at t_0
    _,_,t_0 = FindAngle(landing)
    print("Initial angle            : ",FindAngle(landing)[1])
    R_0,_,v = landing.orient()
    print("Initial angular velocity : ", np.linalg.norm(v)/np.linalg.norm(R_0))

    d_t = 1400#(60*60*24)*0.01

    """
    The following code was part of the process of finding our landing coordinates
    It generates 10 images spaced around the planet (not neccesarily a full orbit).
    It's been commented out so we don't generate 11 image files when running this program
    """
    # coords_list = []
    # N = 10
    # phis = []
    # for n in range(N):

    #     # Finding the position at the current time
    #     r = mission.system.radii[1] * 1000
    #     _,phi,t = FindAngle(landing)

    #     # Defining our coordinate vector
    #     coord_vector = np.array((r,phi,0))

    #     # Finding the vector at time = 0
    #     theta_0 = AngleAtZero(coord_vector,t,p_phi)
    #     coord_vector_0 = np.array((r,theta_0,0))

    #     # saving data
    #     coord_info = (coord_vector,coord_vector_0,t)
    #     coords_list.append(coord_info)

    #     # Taking photo
    #     landing.take_picture(f"Preparation_image_at_{int(t)}_time.xml")

    #     # Updating position
    #     landing.fall(d_t)

    #     phis.append(phi)

    # Simply dropping the thing
    boost = np.array([0,0,0])
    landing.boost(boost)
    landing.fall(d_t+100)
    landing.take_picture(f"Target.xml")

    r           = mission.system.radii[1] * 1000
    R_vec,phi,t = FindAngle(landing)
    # This theta value was found by
    theta       = np.arccos((-1000 * (t-t_0))/np.linalg.norm(R_vec))
    theta       = np.arccos((-1000 * (t-t_0))/np.linalg.norm(R_vec))

    r_0     = r
    phi_0   = AngleAtZero(np.array((r,phi,theta)),t-t_0,p_phi)
    theta_0 = theta

    """ Landing info """
    print("\nOur landing position at time t = 0:")
    print(f"[{r_0},{phi_0},{theta_0}]")
    print(f"Our landing position at time t = {t-t_0}:")
    print(f"[{r},{phi},{theta}]")

    """ Testing info """
    print(f"\n Reinsertion of dt = {t - t_0}: ")
    r,phi,theta = CoordinateAtTime(np.array((r_0,phi_0,theta_0)),t - t_0,p_phi)
    print(f"[{r},{phi},{theta}]")
    P = mission.system.rotational_periods[1]
    print(f"\n Skipping a rotation:")
    r,phi,theta = CoordinateAtTime(np.array((r_0,phi_0,theta_0)),t - t_0 + P,p_phi)
    print(f"[{r},{phi},{theta}]")