# Regular imports
import numpy             as np
import matplotlib.pyplot as plt
from numba import njit
import pickle as pkl 

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

"""
Our seed is [29578].
"""

k = const.k_B * 10 ** (18)

@njit
def ChiSquare(fluxes, noises,l_0,m, T, F_min):

    num_wavelengths = len(fluxes)

    O = np.zeros(3)
    X_min = 10e20

    # Looping over all parameters
    for i_l in range(len(fluxes)):
        for i_t in range(len(T)):
            for i_f in range(len(F_min)):
                
                l     = fluxes[i_l][0]
                t     = T[i_t]
                f_min = F_min[i_f]

                X = 0

                # Looping over all fluxes
                for i in range(len(fluxes)):
                    flux  = fluxes[i][1]
                    noise = noises[i][1]

                    X += (flux - (1 + (f_min - 1)*np.exp(-0.5 * ((3*10**8)* (l - l_0))/(l_0 *((t*k/m)**(1/2))))))/noise

                if X < X_min:
                    X_min = X
                    O[0] = l
                    O[1] = t
                    O[2] = f_min

    return O

njit
def ChiRange(fluxes, noises, L_0,M):

    T = np.linspace(1,1000,1000)
    F_min = np.linspace(0,1,1000)

    O = np.zeros((len(L_0),3))

    for i in range(len(L_0)):

        i_l_min,i_l_max = FindLambda(L_0[0]) # Note that these are indexes

        flux_range = fluxes[i_l_min:i_l_max]
        noise_range = noises[i_l_min:i_l_max]

        print(f"Running optimization for lambda_0 : {L_0[i]}")
        var = ChiSquare(flux_range,noise_range,L_0[i],M[i],T,F_min)
        O[i] = var

    return O

@njit
def FindLambda(l_0):
    d_l = (1/3) * l_0 * 10**(-3)

    l_max = l_0 + d_l
    l_min = l_0 - d_l

    i_max = l_max // 0.0005
    i_min = l_min // 0.0005

    return int(i_min),int(i_max)

if __name__ == "__main__":

    # Ast init
    seed = utils.get_seed('bmthune')
    mission = SpaceMission(seed)  
 
    # Files we will extract data from.
    noise_path = "sigma_noise.txt"
    flux_path = "spectrum_seed78_600nm_3000nm.txt"

    """
    This way of reading the data is very very slow
    so maybe we change it ?
    """

    noises = []
    with open(noise_path,"r") as file:
        for line in file:
            line = file.readline().strip(" ")
            words = line.split("     ")

            wavelength = float(words[0])
            noise = float(words[1])

            noises.append((wavelength,noise))

    fluxes = []
    with open(flux_path,"r") as file:
        for line in file:
            line = file.readline().strip(" ")
            words = line.split("     ")

            wavelength = float(words[0])
            flux = float(words[1])

            fluxes.append((wavelength,flux))

    L_0 = np.array([ 
        632,
        690,
        760,
        720,
        820,
        940,
        1400,
        1600,
        1660,
        2200,
        2340,
        2870
        ])
    
    M   = np.array([ 
        32,
        32,
        32,
        18,
        18,
        18,
        44,
        44,
        16,
        16,
        28,
        44
        ]) * const.m_p
    
    print("Running Optimization")
    O = ChiRange(fluxes,noises,L_0,M)

    print("")

    for i in range(len(O)):
        print(f"Optimization for lambda_0 = {L_0[i]} : lambda = {O[i][0]:20}, T = {O[i][1]:20}, F_min = {O[i][2]:20}")
