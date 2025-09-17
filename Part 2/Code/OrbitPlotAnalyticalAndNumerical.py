# Regular imports
import numpy             as np
import matplotlib.pyplot as plt

# Orbit imports
from OrbitPlotAnalytical import AnalyticalOrbit
from OrbitPlotNumerical import NumericalOrbitFunction

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils

from ast2000tools.space_mission import SpaceMission



# Initializing AST
Seed = utils.get_seed('bmthune')
mission = SpaceMission(Seed)
system = mission.system

# Instantiating Analytical orbit class (and running loop)
OrbitA = AnalyticalOrbit(SemiMajors = system.semi_major_axes,Eccentricities = system.eccentricities,AphelionAngles=system.aphelion_angles)
r_A = OrbitA.Loop()

# Instantiating the Numerical Orbit class (and running the loop)
OrbitN = NumericalOrbitFunction("NumericalOrbitData.npz")

r_N = OrbitN.range(0,OrbitN.TotalTime)

# Initializing plotting
fig, ax = plt.subplots()

# Plotting Orbit data (Analytical)
colors = OrbitA.GetColors()
for i in range(OrbitA.NumPlanets):
    ax.plot(r_A[i][0],r_A[i][1],color = OrbitA.primary)

# Plotting the orbits of the planets
colors = OrbitN.GetColors()
for i in range(len(r_N[0])):
    ax.plot(r_N[0][i],r_N[1][i], color = OrbitN.primary)

# Adding the star (Not to scale)
star = plt.Circle((0, 0), 0.75, color = 'gold')
ax.add_patch(star)

# Adding title, and axis labels
plt.title("Plott av planetenes baner, beregnet \n analytisk (Rød) og numerisk (Lilla)")
plt.xlabel("Posisjon langs x-aksen [AU]")
plt.ylabel("Posisjon langs y-aksen [AU]")

# Making axes equal, and showing plot
plt.axis('equal')
plt.show()