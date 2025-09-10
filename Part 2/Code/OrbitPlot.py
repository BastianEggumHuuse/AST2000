
# AST imports
import ast2000tools.constants as const
import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem

seed = utils.get_seed('bmthune')
system = SolarSystem(seed)

print(system.eccentricities[0])
print(system.semi_major_axes[0])
print(system.semi)

