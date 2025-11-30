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
    phi_0 = coordinate_at_zero[1]
    theta_0 = coordinate_at_zero[2]

    r = r_0
    phi = phi_0 + p_theta * elapsed_time

    coordinate_at_time = np.array((r,phi,theta_0))
    return coordinate_at_time

def AngleAtZero(coordinate_at_time,elapsed_time,p_phi):

    """
    Function that computes the position of a point on the planet at time = 0

    parameters:
    coordinate_at_time : array(float) | coordinates (spherical) of the point we wish to find, at given time
    elapsed_time       : float        | the time which has elapsed
    p_phi            : float        | angular velocity of points on the planet (how fast the planet spins)

    returns:
    float | the angle phi of the point given with coordinate_at_time, but at time = 0
    """

    phi = coordinate_at_time[1]
    
    phi_0 = phi - p_phi * elapsed_time

    return phi_0

def FindAngle(landing_sequence):

    t,r,v = landing_sequence.orient()

    phi = ComputeAngle(r)

    return r,phi,t

def ComputeAngle(r):

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

    return phi

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
    print("Planet angular velocity  : ", p_phi)

    # Finding the time, position, and velocity at t_0
    _,_,t_0 = FindAngle(landing)
    print("Initial angle            : ",FindAngle(landing)[1])
    R_0,_,v = landing.orient()
    print("Initial angular velocity : ", np.linalg.norm(v)/np.linalg.norm(R_0), "\n")

    d_t = 1400#(60*60*24)*0.01

    """
    The following code was part of the process of finding our landing coordinates
    It generates 10 images spaced around the planet (not neccesarily a full orbit).
    It's been commented out so we don't generate 11 image files when running this program
    """
    # coords_list = []
    # phis = []
    # N = 10
    # for n in range(N):

    #   # Finding the position at the current time
    #     r = mission.system.radii[1] * 1000
    #     r_vec_1, phi,t = FindAngle(landing)

    #    # Defining our coordinate vector
    #     coord_vector = np.array((r,phi,0))

    #      # Finding the vector at time = 0
    #     phi_0 = AngleAtZero(coord_vector,t,p_phi)
    #     coord_vector_0 = np.array((r,phi_0,0))

    #     # saving data
    #     coord_info = (coord_vector,coord_vector_0,t)
    #     coords_list.append(coord_info)

    #     # Taking photo
    #     landing.take_picture(f"Preparation_image_at_{int(t)}_time.xml")

    #     # Updating position
    #     landing.fall(d_t)

    #     phis.append(phi)
   
    boost = np.array([0,0,0])
    landing.boost(boost)
    landing.fall(d_t)
    landing.take_picture(f"target_just_before.xml")
    landing.fall(30)
    landing.take_picture(f"target_just_right.xml")

    r           = mission.system.radii[1] * 1000
    r_vec,phi,t = FindAngle(landing)
    theta       =  np.arccos(r_vec[2]/np.linalg.norm(r_vec))

    r_0     = r
    phi_0   = AngleAtZero(np.array((r,phi,theta)),t-t_0,p_phi)
    theta_0 = theta
    
    print("\nOur landing position at time t = 0:")
    print(f"[{r_0},{phi_0},{(theta_0)}]")
    print(f"Our landing position (angles) at time t = {t-t_0}:")
    print(f"[{r},{phi},{theta}]")

    """ Testing info """
    print(f"\nReinsertion of dt = {t - t_0}: ")
    r,phi,theta = CoordinateAtTime(np.array((r_0,phi_0,theta_0)),t - t_0,p_phi)
    print(f"[{r},{phi},{theta}]")
    P = mission.system.rotational_periods[1] * 24 * 60 * 60
    print("\nRotational time : ", P)
    print(f"\nSkipping a rotation:")
    r,phi,theta = CoordinateAtTime(np.array((r_0,phi_0,theta_0)),t - t_0 + P,p_phi)
    print(f"[{r},{phi},{theta}]")
    print(f"[{r},{phi - 2*np.pi},{theta}]")

"""
Output:

Planet angular velocity  :  4.7248043313249e-05
Initial angle            :  3.4578945714340987
Initial angular velocity :  0.05678614419623254

XML file target_just_before.xml was saved in XMLs/.
It can be viewed in MCAst.
XML file target_just_right.xml was saved in XMLs/.
It can be viewed in MCAst.

Our landing position at time t = 0:
[2304594.3015970597,4.794607172841579,1.6083681266857124]
Our landing position (angles) at time t = 1430.0:
[2304594.3015970597,4.862171874779525,1.6083681266857124]

Reinsertion of dt = 1430.0:
[2304594.3015970597,4.862171874779525,1.6083681266857124]

Rotational time :  132982.97382439315

Skipping a rotation:
[2304594.3015970597,11.145357181959112,1.6083681266857124]
[2304594.3015970597,4.862171874779525,1.6083681266857124]
"""