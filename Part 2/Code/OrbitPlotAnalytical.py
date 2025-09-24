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

    def __init__(self,SemiMajors,Eccentricities,AphelionAngles):
        
        # Storing Semi-major axes and Eccentricities for all planets
        self.a = SemiMajors
        self.e = Eccentricities
        self.AphAngles = AphelionAngles
        # Getting number of planets, for looping later
        self.NumPlanets = len(self.a)

        # We wish to go through a full rotation for each planet, so we want an r value for all angles between 0 and 2pi
        # We have chosen 1000 as our number of angles, which gives more than sufficient accuracy
        self.NSteps = 1000
        self.Angles = np.linspace(0,2*np.pi, self.NSteps)
        # Out r-array has the dimention (number of planets, 2, number of angles). This lets us operate on all angles at the same time through vectorization
        self.r = np.zeros((self.NumPlanets,2, len(self.Angles)))

        # Colors we display the different planets with
        self.colors = [[1,0.8,0], [1,0.7,0], [1,0.6,0], [1,0.5,0], [1,0.4,0], [1,0.2,0], [1,0.0,0]]
        # If we want to display the planets with just one color, we use this one instead
        self.primary = [1,0,0]
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
            # COMMENT THAT THIS WAS INCORRECT!!!!!!!!!!!
            r = ((self.a[i]*(1-self.e[i]**2))/(1+self.e[i]*np.cos(self.Angles + np.pi)))

            # Setting first all the x-axis positions, and then all the y-axis-positions.
            # COMMENT THAT THIS WAS INCORRECT!!!!!!!!!!!
            self.r[i][0] = r * np.cos(self.Angles + self.AphAngles[i])
            self.r[i][1] = r * np.sin(self.Angles + self.AphAngles[i])
        # Returning all our positions
        return(self.r)

if __name__ == "__main__":
      
    # Initializing ast2000tools
    seed = utils.get_seed('bmthune')
    mission = SpaceMission(seed)
    system = mission.system

    # Instantiating Analytical orbit class (and running loop)
    Orbit = AnalyticalOrbit(SemiMajors = system.semi_major_axes,Eccentricities = system.eccentricities,AphelionAngles=system.aphelion_angles)
    r = Orbit.Loop()

    # Initializing plotting
    fig, ax = plt.subplots()

    # Plotting Orbit data
    colors = Orbit.GetColors()
    for i in range(Orbit.NumPlanets):
        ax.plot(r[i][0],r[i][1],color = colors[i])

    # Adding the star (not to scale)
    star = plt.Circle((0, 0), 0.75, color = 'gold')
    ax.add_patch(star)

    # Adding title, and axis labels
    plt.title("Plott av planetenes baner, beregnet analytisk")
    plt.xlabel("Posisjon langs x-aksen [AU]")
    plt.ylabel("Posisjon langs y-aksen [AU]")

    # Making axes equal, and showing plot
    plt.axis('equal')
    plt.show()



    # TESTING ZONE!!

    # Defining an epsilon as an upper limit to how inaccurate our tests can be
    eps = 1e-3
    # Also keeping track of how many tests we perform
    TestCount = 0

    '''
    There are three tests we want to perform:
    1) The way we've designed the system, each planet begins in their Aphelion position. We test if this is true, by checking if the initial position has the largest radius.
    2) Similarily, we test if the halwaypoint (time step 500) has the lowest radius (this would be the Perihelion position)
    3) Lastly, we check if the sum of these radii is equal to 2*a, where a is the planet's semi-major axis (retrieved from ast2000tools).
    These two values should be equal, since the sum of the aphelion point and the perihelion point should equal the entire length of the ellipse (along it's longest axis)
    '''

    for i in range(Orbit.NumPlanets):
        
        # Retrieving the positions for time step 0 (aphelion) and 500 (perihelion)
        # Remember that the indexes are [planet][axis][step]
        r_0 = np.linalg.norm(np.array([r[i][0][0],r[i][1][0]]))
        r_pi = np.linalg.norm(np.array([r[i][0][500],r[i][1][500]]))

        # Retrieving the maximum and minimum radii from the orbits
        norms = np.linalg.norm(np.array([r[i][0],r[i][1]]).T,axis = 1, keepdims = True)
        r_min = min(norms)[0]
        r_max = max(norms)[0]

        # Performing tests

        print(f"Planet : {i + 1}")
        # Test 1
        if not (abs(r_0 - r_max) / abs(r_max) < eps):
            raise ValueError(f"Aphelion radius for planet {i+1} not correct!\nValue is {r_0} AU, but should be {r_max} AU")
        print(f"Aphelion radius : {r_max}, First radius : {r_0}, Relative error {abs(r_0 - r_max) / abs(r_max)}")
        # Test 2
        if not (abs(r_pi - r_min) / abs(r_min) < eps):
            raise ValueError(f"Perihelion radius for planet {i+1} not correct!\nValue is {r_pi} AU, but should be {r_min} AU")
        print(f"Perihelion radius : {r_min}, Halfway radius : {r_pi}, Relative error {abs(r_pi - r_min) / abs(r_min) }")
        # Test 3
        if not (abs((r_0 + r_pi) - 2*Orbit.a[i]) / abs(2*Orbit.a[i]) < eps):
            raise ValueError(f"Total ellipse radius for planet {i+1} not correct!\nValue is {r_0 + r_pi} AU, but should be {2*Orbit.a[i]} AU")
        print(f"Sum of Aphelion and Perihelion radii : {r_0 + r_pi}, 2a: {2*Orbit.a[i]}, Relative error {abs(r_0 + r_pi - 2*Orbit.a[i]) / abs(2*Orbit.a[i]) }")
        print("")
        TestCount += 3

    # All is good :)
    print(f"All {TestCount} tests performed, none failed.")
