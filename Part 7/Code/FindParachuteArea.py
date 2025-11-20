# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune

# Regular imports
import numpy             as np
import matplotlib.pyplot as plt
from numba import njit
import pickle as pkl 

from HabitableZones import HabitableZones

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

#Unpickling launch
with open("Mission.pkl", 'rb') as file:
    mission = pkl.load(file)

# Constants from Analytiskeløsninger.py
Zones  = HabitableZones(mission)
Temp_s = Zones.Loop()[1]
T_s    = np.mean(Temp_s)
mu     = 16
r_s    = mission.system.radii[1]*1000
rho_s  = mission.system.atmospheric_densities[1]
P_s    = rho_s * const.k_B * T_s/(mu*const.m_p)
gamma  = 1.4
g      = const.G * mission.system.masses[1]*const.m_sun /(r_s**2)

# Function taken directly from Analytiskeløsninger.py
def rho(r):
    a = P_s**(1-gamma) * T_s**gamma
    C = 5/2 * rho_s**(2/5) +  a*g/gamma * (mu*const.m_p/const.k_B)**gamma * r_s
    rho = (2/5*(-a*g/gamma * (mu*const.m_p/const.k_B)**gamma * r + C))**(5/2)
    return rho

# Calculating velocity for a safe landing
surface_radius   = mission.system.radii[1] * 1000
planet_mass      = mission.system.masses[1] * const.m_sun
lander_mass      = mission.lander_mass
surface_density  = rho(surface_radius) 
planet_rotation  = (2*np.pi) / (mission.system.rotational_periods[1] * 24 * 60 * 60)

radial_velocity  = 3
air_velocity     = surface_radius * planet_rotation
drag_velocity_2  = radial_velocity**2# + air_velocity**2

area = (2 * const.G * planet_mass * lander_mass) / (rho_s * drag_velocity_2 * surface_radius**2)
print(f"Neccesary area: {area} m^2")

