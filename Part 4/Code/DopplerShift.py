# Ikke brukt kodemal!!
# Skrevet av Bendik Thune og Bastian Eggum Huuse

import sys
import numpy as np
from numba import njit
import pickle as pkl 


import ast2000tools.constants as const
import ast2000tools.utils     as utils
import GeneralizedLaunch as Lan 
from ast2000tools.space_mission import SpaceMission

@njit
def Doppler2V(dLambda, Lambda_0):
    """
    Finds the radiall speed of one star by messuring the dopple effect
    Inputs:
    dLambda: Messured Doppler shift (nm)
    Lambda_0: Original spectral line(nm)
    Outputt: Velocity in radial direction gotten from dopplershift (m/s)
    """
    v_r = dLambda/Lambda_0 * const.c #the formula for dopler shift
    return v_r


@njit
def RelVel(dLambda_1, dLambda_2,Lambda_0):
    """
    Creates a tupple with the velocity relativ to two stars
    Inputs:
    dLambda_1: Messured dopplereffect of star 1 (nm)
    dLambda_2: Messured dopplereffect of star 2 (nm)
    Lambda_0: Original spectral line (nm)
    Outputt: Vellocity given in the basis of star 1 and 2  (m/s)
    """

    v_1 = Doppler2V(dLambda_1,Lambda_0) 
    v_2 = Doppler2V(dLambda_2, Lambda_0)
    return (v_1,v_2)

@njit
def BasisShift(phi_1,phi_2, v):
    """
    Changes basis from star basis to xy
    Inputs:
    phi_1: angle to star 1
    phi_2: angle to star 2
    v: velocity given in basis of speed away from star 1 and 2
    Outputt: velocity in xy basis
    """
    #Here we change the basis from stars to x,y notice adding pi, this is because the angles gives us the direction toward the stars
    #while the basis goes the other way
    v_x = np.cos(phi_1 + np.pi)*v[0] + np.cos(phi_2+np.pi) *v[1]
    v_y = np.sin(phi_1 + np.pi)*v[0] + np.sin(phi_2 +np.pi) *v[1]
    return v_x, v_y

@njit
def RockVel(v_r,v_s):
    """
    Find the relativ velocity of rocket to sun
    Inputs:
    v_r: vellocity of rocket compared to stars 
    v_s: vellocity of sun compared to stars (v_s and v_r must have same unit)
    outputt: Vellocity of rocket compared to sun (same unit as sent in)
    """
    v_x = v_r[0] - v_s[0]
    v_y = v_r[1] - v_s[1]
    return v_x, v_y #not nesesarly in xy basis

@njit
def Main(dLambda, Lambda_0, phi_1, phi_2):
    """
    Runns all functions above to find the vellocity 
    Inputs:
    dLambda: an array of length 4 that has the dopplershifts from each star on the sun and rocket (in nm)
    Lambda_0: Original spectral line in nm
    phi_1: angle to star 1
    phi_2: angle to star 2
    Outputt: Velocity of rocket in xy basis compared to sun (AU/Y)
    """
    #first i find the velocity of rocket at sun in star basis
    v_s = RelVel(dLambda[0],dLambda[1],Lambda_0)
    v_r = RelVel(dLambda[2],dLambda[3],Lambda_0)
    
    #Find the velocity of rocket compared to sun in star basis
    V = RockVel(v_r, v_s)
    
    #change the basis to xy basis
    V_r = BasisShift(phi_1,phi_2, V)
    
    #returns the velocity converted to AU/Y
    return V_r[0]*60*60*24*365/const.AU, V_r[1]*60*60*24*365/const.AU

def tests(mission):
    """
    Test function that runns two tests, if vel of rocet is equal to sun and if equal to stars compared to sun, 
    returns nothing if sucsesfull
    """
    epsilon = 1e-3 #A small epsilon to check if a value is aprox zero
    lambda_0 = 656.3 #wavelength of spectral line of H_alpha, goten from the problem sett
    
    #Sets the lambdas where rocket has 0 vel compared to sun, therfore will have same dopplerdisplacement as the sun
    lambda_1, lambda_2 = mission.star_doppler_shifts_at_sun 
    dlambda = (lambda_1, lambda_2,lambda_1, lambda_2)
    #gets angles
    phi_1, phi_2 = mission.star_direction_angles
    phi_1, phi_2 = (np.deg2rad(phi_1), np.deg2rad(phi_2))
    #finds rocket speed should be zero
    V_r = Main(dlambda, lambda_0, phi_1, phi_2)
    #Raises error if not zero
    if not(abs(V_r[0]) < epsilon and abs(V_r[1]) < epsilon):
        
        raise ValueError(f'Expected V = 0 got V = {V_r}AU/Y')
    
    #sets the shifts in the rocket to zero so it has same speed as stars
    dlambda = (lambda_1, lambda_2, 0, 0)
    
    #should then get V_r equal to negative v_s
    V_r = Main(dlambda, lambda_0, phi_1, phi_2)
    v_s = BasisShift(phi_1,phi_2, RelVel(dlambda[0],dlambda[1],lambda_0))
    v_s = np.array(v_s) *60*60*24*365/const.AU

    #raises error if V_r not -v_s
    if not((abs(V_r[0] + v_s[0])) < epsilon and abs(V_r[1] + v_s[1]) < epsilon):
        raise ValueError(f'Expected V = -V_sun, got V+v_sun = {V_r[0] + v_s[0],V_r[1] + v_s[1]} AU/Y')
    

 
if __name__ == '__main__':
    #gets the mission that has been saved using pickle 
    #THis saves having to simulate the launch every time
    with open("Mission.pkl", 'rb') as file:
        mission = pkl.load(file)

    lambda_0 = 656.3 #Hydrogen spectral line
    lambda_1, lambda_2 = mission.star_doppler_shifts_at_sun #gets dopler form sun
    lambda_3, lambda_4 = mission.measure_star_doppler_shifts() #gets doper from rocket after launch
    dlambda = (lambda_1,lambda_2,lambda_3,lambda_4) #creates the array with all doppler shifts
  

    phi_1, phi_2 = mission.star_direction_angles #gets angles
    phi_1, phi_2 = (np.deg2rad(phi_1), np.deg2rad(phi_2)) #converts them to radians, why does the mission give them in degrees???

    tests(mission)#runs the test function
    
    V_r = Main(dlambda, lambda_0, phi_1, phi_2) #runs main
    print(V_r) #returns the vellocity after launch

    """
    Runtime example:
    (-6.367980083713005, -0.692754340951314)

    """