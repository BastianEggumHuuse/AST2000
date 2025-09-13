# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune


# Regular imports
import numpy             as np
import matplotlib.pyplot as plt


# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils

from ast2000tools.space_mission import SpaceMission

class AnalyticalOrbit:

    def __init__(self,SemiMajors,Eccentricities):
        
        # Storing Semi-major axes and Eccentricities for all planets
        self.a = SemiMajors
        self.e = Eccentricities
        # Getting number of planets, for looping later
        self.NumPlanets = len(self.a)

        # We wish to go through a full rotation for each planet, so we want an r value for all angles between 0 and 2pi
        # We have chosen 1000 as our number of angles, which gives more than sufficient accuracy
        self.Angles = np.linspace(0,2*np.pi, 1000)
        # Out r-array has the dimention (number of planets, 2, number of angles). This lets us operate on all angles at the same time through vectorization
        self.r = np.zeros((self.NumPlanets,2, len(self.Angles)))

        # Colors we display the different planets with
        self.colors = [[1,0.8,0], [1,0.7,0], [1,0.6,0], [1,0.5,0], [1,0.4,0], [1,0.2,0], [1,0.0,0]]
        # If we want to display the planets with just one color, we use this one instead
        self.primary = [1,0.2,0]
        #self.colors = ["red", 'darkorange', 'mediumblue', 'limegreen', 'purple', 'darkviolet', 'gold']

    def GetColors(self):

        """Method that returns the colors of the planet orbits in the correct order."""

        # Color index | Star index
        # 0           | 0
        # 1           | 1
        # 2           | 4
        # 3           | 3
        # 4           | 5
        # 5           | 6
        # 6           | 2

        # The planets aren't listed by distance from the sun, so we have to move the colors around a bit so that they match.
        # See above table for detailed indexes.
        return([self.colors[0],self.colors[1],self.colors[4],self.colors[3],self.colors[5],self.colors[6],self.colors[2]])

    def Loop(self):
        
        """Looping through the planets, and calculating their positions for one entire rotation"""

        for i in range(self.NumPlanets):

            # We find the distance from the star at for all the angles (vectorized)
            r = ((self.a[i]*(1-self.e[i]**2))/(1+self.e[i]*np.cos(self.Angles)))

            # Setting first all the x-axis positions, and then all the y-axis-positions.
            self.r[i][0] = r * np.cos(self.Angles)
            self.r[i][1] = r * np.sin(self.Angles)

        # Returning all our positions
        return(self.r)

if __name__ == "__main__":
      
    # Initializing ast2000tools
    seed = utils.get_seed('bmthune')
    mission = SpaceMission(seed)
    system = mission.system

    # Instantiating Analytical orbit class (and running loop)
    Orbit = AnalyticalOrbit(SemiMajors = system.semi_major_axes,Eccentricities = system.eccentricities)
    r = Orbit.Loop()

    # Initializing plotting
    fig, ax = plt.subplots()

    # Plotting Orbit data
    colors = Orbit.GetColors()
    for i in range(Orbit.NumPlanets):
        ax.plot(r[i][0],r[i][1],color = colors[i])

    # Adding the star
    star = plt.Circle((0, 0), 1, color = 'gold')
    ax.add_patch(star)

    # Adding title, and axis labels
    plt.title("Plott av planetenes baner, beregnet analytisk")
    plt.xlabel("Posisjon langs x-aksen [AU]")
    plt.ylabel("Posisjon langs y-aksen [AU]")

    # Making axes equal, and showing plot
    plt.axis('equal')
    plt.show()