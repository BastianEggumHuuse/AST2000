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

# Initializing ast2000tools
seed = utils.get_seed('bmthune')
mission = SpaceMission(seed)
system = mission.system

Func = NumericalOrbitFunction("NumericalOrbitData.npz")

# Instantiating Analytical orbit class (and running loop)
Orbit = AnalyticalOrbit(SemiMajors = system.semi_major_axes,Eccentricities = system.eccentricities,AphelionAngles=system.aphelion_angles)
r = Orbit.Loop()

# Initializing plotting
fig, ax = plt.subplots()

# Plotting Orbit data
colors = Orbit.GetColors()
for i in range(Orbit.NumPlanets):
    ax.plot(r[i][0],r[i][1],color = colors[i])

r = Func.range(0,Func.TotalTime)

# Plotting the orbits of the planets
colors = Func.GetColors()
for i in range(len(r[0])):

    r_time = (2*np.pi)*((system.semi_major_axes[i]**3)/(const.G_sol*system.star_mass + system.masses[i]))**(1/2)
    r_t = Func.range(0, r_time)
    # Finding radius at aphelion and perihelion COMMENT THIS!!!!!!!!!
    t_aph = 0
    t_per = 0
    r_aph = np.linalg.norm(np.array([r_t[0][i][0],r_t[1][i][0]]))
    r_per = r_aph

    for j in range(len(r_t[0][0])):

        r_a = np.linalg.norm(np.array([r_t[0][i][j],r_t[1][i][j]]))

        if(r_a < r_per):
            r_per = r_a
            t_per = j * Func.dt
        if(r_a > r_aph):
            r_aph = r_a
            t_aph = j * Func.dt

    if(t_per < t_aph):
        t_per = r_time + t_per

    print(t_aph,t_per)

    r_c = Func.range(t_aph, t_per)

    ax.plot(r_c[0][i],r_c[1][i], color = colors[i])

print(system.aphelion_angles)

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