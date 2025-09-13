# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune


# Regular imports
import numpy             as np
import matplotlib.pyplot as plt


# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils

from ast2000tools.space_mission import SpaceMission

class NumericalOrbit:
    def __init__ (self,mission, const, TotalTime, StepsPerYear, InitialPos, InitialVel):
        
        # Storing mission, system, and constants
        self.mission = mission
        self.system = mission.system 
        self.const = const

        # Storing gravitational constant and mass of out star
        self.G = const.G_sol
        self.SM = self.system.star_mass

        # Storing TotalTime and StepsPerYear (Time is measured in years)
        self.T = TotalTime
        self.YSteps = StepsPerYear

        # Setting deltatime and N timesteps
        self.dt = 1/StepsPerYear
        self.NSteps = int(TotalTime/self.dt)

        # Setting initial positions and velocities. These are input with dimentions (2,Num_planets), meaning the first dimention contains two elements, 
        # which each contain all x-values, and all y-values.
        # Since we want to be able to operate on all the planets simultaneously, we want the dimentions (7,Num_planets), 
        # meaning the first dimention contains Num_planets elements, which contain a single x and y element.
        self.r0 = InitialPos.T
        self.v0 = InitialVel.T
        
        # Creating our arrays. Again, we want the dimentions (Num_planets,2) to operate on all planets at the same time.
        # Since we now want to go through time as well, we add this dimention at the beginning, giving us the dimentions (N_timesteps,Num_planets,2)
        self.NumPlanets = len(self.r0)
        self.r = np.zeros((self.NSteps, self.NumPlanets, 2))
        self.v = np.copy(self.r)
        self.a = np.copy(self.r)

        # We create a time array with linearly spaced values between 0 and T, with NSteps elements.
        self.t = np.linspace(0,self.T,self.NSteps) 

        # Setting the initial values for the arrays. r0 and v0 are already defined, but accelerations has to be calculated.
        # We use Newtons law of gravitation to get the accelerations.
        self.r[0] = self.r0
        self.v[0] = self.v0
        self.a[0] = ((-self.G*self.SM)/self.norm(self.r[0])**2)*self.hat(self.r[0])

        # Colors we display the different planets with
        self.colors = [[0,0,1], [0.3,0,1], [0.4,0,1], [0.5,0,1], [0.6,0,1], [0.7,0,1], [0.8,0,1]]
        # If we want to display the planets with just one color, we use this one instead
        self.primary = [0.5,0,1]

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

    def norm(self, v):

        """
        Support method that gets the norms of an array of vectors.
        Used by us mainly to get |r| for all the planets at once.
        """

        return np.linalg.norm(v, axis=1, keepdims=True)
    
    def hat(self, v):
        
        """
        Support method that returns the unit vector for a given vector
        Used by us mainly to get r_hat for the gravitational acceleration.
        (Compatible with vectorized code)
        """

        return v/self.norm(v)

    def timestep(self,i):
        
        """
        Performs one time step within the orbit simulation. 
        The parameter [i] refers to which timestep we are currently at,
        meaning that 0 <= i < NSteps - 1
        """

        #self.r[i] = self.r[i-1] + self.v[i-1]*self.dt + 0.5 * self.a[i-1]*self.dt**2 
        #self.a[i] = ((-self.G*self.SM)/self.norm(self.r[i])**2)*self.hat(self.r[i])
        #self.v[i] = self.v[i-1] + 0.5*(self.a[i-1] + self.a[i]) *self.dt

        # Performing leapfrog integration. Note that because of our vectorization, every operation that happens here happens for all planets at once.
        self.r[i+1] = self.r[i] + self.v[i]*self.dt + 0.5 * self.a[i]*self.dt**2 
        self.a[i+1] = ((-self.G*self.SM)/self.norm(self.r[i+1])**2)*self.hat(self.r[i+1])
        self.v[i+1] = self.v[i] + 0.5*(self.a[i] + self.a[i+1]) *self.dt
    
    def loop(self):

        """
        Performs the entire loop of the orbit simulation. Since the timestep refers to i+1, we want to skip the final step, to avoid an indexing error.
        Therefore, we loop from 0 (included) to (NSteps - 1) (Excluded).
        The method also returns the arrays r, v, a, and t.
        """

        # Loop
        for i in range(0,self.NSteps-1):
            # Performing the step
            self.timestep(i)     

        # Returning our arrays for plotting. We transpose the arrays again before returning,
        # which gives us the arrays with dimentions (2,NumPlanets,NSteps), letting us plot for all timesteps.
        return(self.r.T,self.v.T,self.a.T,self.t)
    

if __name__ == "__main__":

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

    # Instantiating the Numerical Orbit class (and running the loop)
    Orbit = NumericalOrbit(mission = mission,const = const, TotalTime = TotalTime, StepsPerYear = 10000, InitialPos = R0, InitialVel = V0)
    r,v,a,t = Orbit.loop()

    # Initializing plotting
    fig, ax = plt.subplots()

    # Plotting the orbits of the planets
    colors = Orbit.GetColors()
    for i in range(Orbit.NumPlanets):
        ax.plot(r[0][i],r[1][i], color = colors[i])

    # Adding the star
    star = plt.Circle((0, 0), 1, color = 'gold')
    ax.add_patch(star)

    # Adding title, and axis labels
    plt.title("Plott av planetenes baner, beregnet numerisk")
    plt.xlabel("Posisjon langs x-aksen [AU]")
    plt.ylabel("Posisjon langs y-aksen [AU]")

    # Making axes equal, and showing plot
    plt.axis('equal')
    plt.show()

    #mission.verify_planet_positions(TotalTime,r)
