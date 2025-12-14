import ast2000tools.utils as utils
from   ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as const
from   ast2000tools.relativity import RelativityExperiments

# Initializing AST 
Seed = utils.get_seed('bmthune')

Experiment = RelativityExperiments(Seed)

Experiment.antimatter_spaceship(1)