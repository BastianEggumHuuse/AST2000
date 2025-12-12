import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import  ast2000tools.constants as const
from ast2000tools.relativity import RelativityExperiments

# Initializing AST 
Seed = utils.get_seed('bmthune')

# Initializing Experiments
Experiments = RelativityExperiments(Seed)

#Experiments.spaceship_duel(1)
#Experiments.cosmic_pingpong(1)
#Experiments.twin_paradox(1)
Experiments.spaceship_race(1)