# Regular imports
import numpy             as np
import matplotlib.pyplot as plt

# Orbit imports
from OrbitPlotNumerical import NumericalOrbitFunction

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils

from ast2000tools.space_mission import SpaceMission

Func = NumericalOrbitFunction("NumericalOrbitData.npz")

r = Func.range(0,Func.TotalTime)

# Initializing plotting
fig, ax = plt.subplots()

# Plotting the orbits of the planets
colors = Func.GetColors()
for i in range(len(r[0])):
    ax.plot(r[0][i],r[1][i], color = colors[i])

# Adding the star (not to scale)
star = plt.Circle((0, 0), 0.75, color = 'gold')
ax.add_patch(star)

# Adding title, and axis labels
plt.title("Plott av planetenes baner, beregnet numerisk")
plt.xlabel("Posisjon langs x-aksen [AU]")
plt.ylabel("Posisjon langs y-aksen [AU]")

# Making axes equal, and showing plot
plt.axis('equal')
plt.show()