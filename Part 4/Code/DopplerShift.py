# Ikke brukt kodemal!!
# Skrevet av Bastian Eggum Huuse

import sys
import numpy as np
from numba import njit

import ast2000tools.constants as const
import ast2000tools.utils     as utils
import GeneralizedLaunch as Lan 
from ast2000tools.space_mission import SpaceMission

@njit
def Doppler2V(dLambda, Lambda_0):
    v_r = dLambda/Lambda_0 * const.c
    return v_r


@njit
def RelVel(dLambda_1, dLambda_2,Lambda_0):
    v_1 = Doppler2V(dLambda_1,Lambda_0) 
    v_2 = Doppler2V(dLambda_2, Lambda_0)
    return (v_1,v_2)

@njit
def BasisShift(phi_1,phi_2, v):
    v_x = np.cos(phi_1 + np.pi)*v[0] + np.cos(phi_2+np.pi) *v[1]
    v_y = np.sin(phi_1 + np.pi)*v[0] + np.sin(phi_2 +np.pi) *v[1]
    return v_x, v_y

@njit
def RockVel(v_r,v_s):
    v_x = v_r[0] - v_s[0]
    v_y = v_r[1] - v_s[1]
    return v_x, v_y

@njit
def Main(dLambda, Lambda_0, phi_1, phi_2):
    v_s = RelVel(dLambda[0],dLambda[1],Lambda_0)
    v_r = RelVel(dLambda[2],dLambda[3],Lambda_0)
    
    V = RockVel(v_r, v_s)
    
    V_r = BasisShift(phi_1,phi_2, V)
    
    
    return V_r[0]*60*60*24*365/const.AU, V_r[1]*60*60*24*365/const.AU
 
if __name__ == '__main__':
    seed = utils.get_seed('bmthune')
    mission = SpaceMission(seed)

    Lan.main(mission, 1)

    lambda_0 = 656.3 #Hydrogen spectral line
    lambda_1, lambda_2 = mission.star_doppler_shifts_at_sun
    lambda_3, lambda_4 = mission.measure_star_doppler_shifts()
    dlambda = (lambda_1,lambda_2,lambda_3,lambda_4)
  

    phi_1, phi_2 = mission.star_direction_angles
    phi_1, phi_2 = (np.deg2rad(phi_1), np.deg2rad(phi_2))
    print(phi_1,phi_2)
    print(mission._velocity_after_launch)
    V_r = Main(dlambda, lambda_0, phi_1, phi_2)
    print(V_r)