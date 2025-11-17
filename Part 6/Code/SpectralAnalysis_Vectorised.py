# Bruker ikke kodemal!!!
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

"""
Our seed is [29578].
"""

k = const.k_B

def ChiSquareOptimization(fluxes, noises,l_0,m,step_length = 1):

    print("-------------------------------------")
    print(f"Running Algorithm for l_0 = {l_0:5} nm")

    # Getting ranges:
    d_l = (1/3) * l_0 * 10**(-3)

    l_max = l_0 + d_l
    l_min = l_0 - d_l

    i_max = int(((l_max - 600) / (2999.9998 - 600)) * len(fluxes))
    i_min = int(((l_min - 600) / (2999.9998 - 600)) * len(fluxes))

    # Getting flux and noise ranges
    flux_range  = fluxes[i_min:i_max]
    noise_range = noises[i_min:i_max]

    print("Observational ranges computed.")

    # Setting range extremes
    sigma_min = (l_0/const.c) * ((const.k_B * 1)/m)**0.5
    sigma_max = (l_0/const.c) * ((const.k_B * 700)/m)**0.5
    f_min_min = 0.5
    f_min_max = 1

    # Getting value range of parameters
    lambda_center_values = flux_range[:,0][::step_length]
    sigma_values         = np.linspace(sigma_min,sigma_max,700)[::step_length]
    f_min_values         = np.linspace(f_min_min,f_min_max,500)[::step_length]
    
    print("Parameter values computed.")

    # Creating grid
    parameter_grid       = np.meshgrid(np.arange(len(lambda_center_values)),np.arange(len(sigma_values)),np.arange(len(f_min_values)))
    lambda_center_range  = parameter_grid[0].flatten()
    sigma_range          = parameter_grid[1].flatten()
    f_min_range          = parameter_grid[2].flatten()

    # Deleting parameter_grid to save memory
    del parameter_grid

    print("Parameter grid created.")

    # Actually optimizing
    X = ChiSquareSum(flux_range,
                     noise_range,
                     lambda_center_values[lambda_center_range],
                     sigma_values[sigma_range],
                     f_min_values[f_min_range]
                     )
    
    print("X sums calculated")
    
    # Finding index of minimum X
    min_X = min(X)
    min_X_index = np.where(X == min_X)[0]
    min_X_index = min_X_index[0]

    # Finding corresponding parameter set
    lambda_center = lambda_center_values[lambda_center_range][min_X_index]
    sigma         = sigma_values[sigma_range][min_X_index]
    temperature   = (m/const.k_B) * ((sigma * const.c)/l_0)**2
    f_min         =f_min_values[f_min_range][min_X_index]

    print(f"Printing results for l_0 = {l_0:5} nm :")

    print("-------------------------------------\n")
    
    print(f"lambda_center : {lambda_center:20} nm")
    print(f"temperature   : {temperature:20} K")
    print(f"f_min         : {f_min:20} flux")

    print("-------------------------------------\n")

def ChiSquareSum(fluxes,noises,lambda_centers,sigmas,f_mins):

    num_fluxes = len(fluxes)
    num_chis   = len(lambda_centers)

    X = np.zeros(num_chis)

    for i in range(num_fluxes):

        lambda_i = fluxes[i][0]
        flux_i   = fluxes[i][1]
        noise_i  = fluxes[i][1]

        X += ((flux_i - (1 + (f_mins - 1)*np.exp(-((lambda_i - lambda_centers)**2)/((2*sigmas**2))))) / noise_i)**2

    return X

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

    fluxes = np.loadtxt(flux_path)
    noises = np.loadtxt(noise_path)

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

    for i in range(len(L_0)):
        ChiSquareOptimization(fluxes,noises,L_0[i],M[i],step_length = 10)

    print("")

    # for i in range(len(O)):
    #     print(f"Optimization for lambda_0 = {L_0[i]:5} : lambda = {O[i][0]:15} nm, delta_lambda = {L_0[i] - O[i][0]:15} nm, T = {O[i][1]:25} K, F_min = {O[i][2]:25}")
