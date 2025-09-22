# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune

# Regular imports
import numpy             as np
import matplotlib.pyplot as plt
from sys import exit

# Orbit imports
from OrbitPlotNumerical import NumericalOrbitFunction

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

class KeplerArea:

    def __init__(self,Filename,PlanetIndex):

        '''
        Class that gets the area swept out by a portion of a planet's orbit.

        Filename    : String | The name of the file storing the numerical orbits.
        PlanetIndex : int    | The index representing the planet we're interested in

        returns     : self
        '''

        # Storing the parameters
        self.Filename = Filename
        self.PlanetIndex = PlanetIndex

        # Instance of NumericalOrbitFunction, which we use to get the numerical positions of the planets.
        self.PositionFunction = NumericalOrbitFunction(self.Filename)

        # Grabbing Variables from PositionFunction (So we don't have to refer to the instance every time)
        self.TotalTime  = self.PositionFunction.TotalTime
        self.dt         = self.PositionFunction.dt
        self.NumSteps   = self.PositionFunction.NumSteps
        self.OrbitTimes = self.PositionFunction.OrbitTimes

    def __call__(self, t0 = 0, t1 = 1):

        """
        Method that runs through and gathers data from NumericalOrbitFunction
        (for the one planet we're interested in)

        t0 : float | Start time in years
        t1 : float | End time in years

        returns : 
        float      | The area swept out by the orbit between t0 and t1
        List       | List of all the position vectors collected from NumericalOrbitFunction
        float      | The total distance traveled between t0 and t1
        float      | The average velocity of the planet within the time interval
        """
        
        # Initializing variables
        R = [] # List of positions
        A = 0  # Area
        S = 0  # Distance
        t = t0 # current time

        while t < t1:
            # Getting the position vector at current t
            r0 = self.PositionFunction(t,self.PlanetIndex)
            # Adding this position to list
            R.append(r0)
            # Updating t
            t += self.dt
            # Getting next position
            r1 = self.PositionFunction(t,self.PlanetIndex)

            #Calculating the distance between each vector and summing them up,
            dS = np.linalg.norm(r1-r0)
            S += dS

            # Finding the area of this small segment of the orbit.
            # Since the segment is so small, it's approximately equal to a triangle
            # We can then use the cross product of the position vectors to get the area!
            # (when using the cross product with vectors of dimention 2,
            # np.cross assumes they are three-dimentional, with z value 0,
            # and then it returns the z-coordinate of the resulting vector,
            # since both x and y will be 0).
            dA = 0.5 * abs(np.cross(r0,r1))
            # Adding to total area
            A += dA

        # Calculating average velocity for this period.
        V = S/(t1- t0)
        
        # Returning all values
        return A, R, S, V
    
        
    
    
if __name__ == "__main__":

    # Initializing AST
    Seed = utils.get_seed('bmthune')
    mission = SpaceMission(Seed)
    system = mission.system

    # We are interested in our own planet, which is the first one
    Planet_index = 0

    # Initializing Area Function
    AreaFunction = KeplerArea("NumericalOrbitData.npz",Planet_index)
    OrbitTime = AreaFunction.OrbitTimes[Planet_index]

    # Finding time at aphelion and perihelion. We know that our planet begins it's orbit in the aphelion position,
    # and therefore that it reaches the perihelion position after half a rotation.
    # We also move start at OrbitTime, instead of 0, such that we never end up with a negative time
    t_aph = OrbitTime + 0
    t_per = OrbitTime + 0.5 * OrbitTime
    # We choose a small interval to check around
    t_interval = 0.05

    # Finding the positions for the entire rotation
    R = AreaFunction.PositionFunction.range(AreaFunction.dt,OrbitTime)   

    # Getting area swept out by orbit, positional vectors, total distance traveled, and average velocity
    A_aph, R_aph, S_aph, V_aph = AreaFunction(t_aph - t_interval,t_aph + t_interval)
    A_per, R_per, S_per, V_per= AreaFunction(t_per - t_interval,t_per + t_interval)

    # Adding the origin to both ends of this array, turning it into a polygon that we can draw.
    Origin = np.zeros((1,2))
    R_aph_polygon = np.concatenate((Origin,R_aph,Origin),0)
    R_per_polygon = np.concatenate((Origin,R_per,Origin),0)
    
    # Printing info (The end of B.1.a, B.1.b, B.1.c)
    print(f"Aphelion area : {A_aph:.5f} AU^2 | Perihelion area : {A_per:.5f} AU^2 | Ratio = {A_aph/A_per:.5f}")
    print(f"Aphelion Distance : {S_aph:.5f} AU | Perihelion Distance : {S_per:.5f} AU")
    print(f"Aphelion Average vel : {V_aph:.5f} AU/Y | Perihelion Average vel : {V_per:.5f} AU/Y ")

    # Initializing plotting
    fig, ax = plt.subplots()

    # Plotting orbits
    ax.plot(R[0][0],R[1][0],color = "green")
    ax.fill(R_aph_polygon[:,0],R_aph_polygon[:,1],color = "royalblue")
    ax.fill(R_per_polygon[:,0],R_per_polygon[:,1],color = "red")
    
    # Making axes equal, and showing plot
    plt.axis('equal')
    plt.show()

    # From here on out we are doing B.2
    # Getting the position function (Gives us position at a time t)
    r = AreaFunction.PositionFunction
    # Making an array with an index for each planet.
    K = np.zeros(7)
    # Defining epsilon and start time
    eps = 1e-3
    t_0 = 1000

    # Looping over all planets
    for p in range(7):
        # Getting the position at the first timestep
        r_0 = r(0,p)

        # We start the simulation a little bit into the orbit, so we don't just grab the starting position :)
        for i in range(t_0,r.NumSteps):
            # Getting the current position
            r_i = r(i*r.dt,p)

            # Checking if relative difference between r_0 and r_i is small enough
            # If it is, we know that we're ca in the same place the orbit began.
            if(np.linalg.norm(r_i - r_0)/np.linalg.norm(r_0) < eps):
                # Calculating the constant
                K[p] = ((i*r.dt)**2)/(system.semi_major_axes[p]**3)
                break

    # Printing info about K
    print(f"\nThe array of constants: {K}")
    print(f"Relative difference between max and min: {abs((max(K) - min(K))/min(K)):.6f}\n")

    # Comparing to the constant from Newtons formula:
    diff = []
    for p in range(7):
        k = K[p]
        # Newtons variant is t^2 = (4pi^2 / G(m1 + m2))*a^3
        # This gives us t^2 / a^3 = (4pi^2 / G(m1 + m2)), where (t^2 / a^3) = K
        Newton_Variant_Constant = (4*np.pi**2)/(const.G_sol*(system.star_mass + system.masses[p]))
        # Getting the relative difference between our value k and the newton_variant_constant
        diff.append(abs((k - Newton_Variant_Constant)/Newton_Variant_Constant))
    
    # Printing the mean of the differences
    print(f"Mean difference between calculated constant and constant aquired from Newtons variant: {np.mean(diff):.7f}")
