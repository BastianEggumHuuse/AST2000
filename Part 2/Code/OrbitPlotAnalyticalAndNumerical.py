# Regular imports
import numpy             as np
import matplotlib.pyplot as plt

# Orbit imports
from OrbitPlotAnalytical import AnalyticalOrbit
from OrbitPlotNumerical import NumericalOrbit

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils

from ast2000tools.space_mission import SpaceMission



# Initializing AST
Seed = utils.get_seed('bmthune')
mission = SpaceMission(Seed)
system = mission.system

# Getting initial conditions
R0 = system.initial_positions
V0 = system.initial_velocities
# Calculating the time the simulation will run. Here we assume that the orbit is a perfect circle, which it isn't, but it's very close.
# To make sure we pass the 20 rotations mark, we multiply the time with 1.1
TotalTime = np.linalg.norm(R0.T[0]) *2*np.pi/np.linalg.norm(V0.T[0])*20*1.1

# Instantiating Analytical orbit class (and running loop)
OrbitA = AnalyticalOrbit(SemiMajors = system.semi_major_axes,Eccentricities = system.eccentricities)
r_A = OrbitA.Loop()

# Instantiating the Numerical Orbit class (and running the loop)
OrbitN = NumericalOrbit(mission = mission,const = const, TotalTime = TotalTime, StepsPerYear = 10000, InitialPos = R0, InitialVel = V0)
r_N,v_N,a_N,t_N = OrbitN.loop()

# Initializing plotting
fig, ax = plt.subplots()

# Plotting Orbit data (Analytical)
colors = OrbitA.GetColors()
for i in range(OrbitA.NumPlanets):
    ax.plot(r_A[i][0],r_A[i][1],color = OrbitA.primary)

# Plotting the orbits of the planets
colors = OrbitN.GetColors()
for i in range(OrbitN.NumPlanets):
    ax.plot(r_N[0][i],r_N[1][i], color = OrbitN.primary)

# Adding the star
star = plt.Circle((0, 0), 1, color = 'gold')
ax.add_patch(star)

# Adding title, and axis labels
plt.title("Plott av planetenes baner, beregnet \n analytisk (Rød) og numerisk (Lilla)")
plt.xlabel("Posisjon langs x-aksen [AU]")
plt.ylabel("Posisjon langs y-aksen [AU]")

# Making axes equal, and showing plot
plt.axis('equal')
plt.show()