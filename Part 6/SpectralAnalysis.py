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

if __name__ == "__main__":

    # Ast init
    seed = utils.get_seed('bmthune')
    mission = SpaceMission(seed)  

    # Files we will extract data from.
    noise_path = "sigma_noise.txt"
    flux_path = "spectrum_seed78_600nm_3000nm.txt"

    print(seed)