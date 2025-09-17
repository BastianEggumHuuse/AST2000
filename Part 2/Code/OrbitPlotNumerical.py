# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune


# Regular imports
import numpy             as np
import matplotlib.pyplot as plt

from OrbitPlotAnalytical import AnalyticalOrbit

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils

from ast2000tools.space_mission import SpaceMission

class NumericalOrbit:
    def __init__ (self,mission, const, TotalTime, StepsPerYear, InitialPos, InitialVel):
        
        """
        Class representing a numerical orbit. The class itself simulates the orbit, when the method loop is called.

        mission      : instance of the SpaceMission class
        const        : instance of the const package
        TotalTime    : float            | the total time the simulation will run for
        StepsPerYear : int              | how many timesteps the simulation will run per year
        InitialPos   : Array(float)     | Array containing initial positions for all planets in the system
        InitialVel   : Array(float)     | Array containing initial positions for all planets in the system

        returns      : self
        """

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

        """
        Method that returns the colors of the planet orbits in the correct order.
        
        returns : list | list of planet colors in correct order
        """

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

        v       : Vector we wish to get the norm of

        returns : float | norm of v
        """

        return np.linalg.norm(v, axis=1, keepdims=True)
    
    def hat(self, v):
        
        """
        Support method that returns the unit vector for a given vector
        Used by us mainly to get r_hat for the gravitational acceleration.
        (Compatible with vectorized code).

        v       : Array(float) | Vector we wish to get the unit vector for.

        returns : Array(float) | unit vector of v
        """

        return v/self.norm(v)

    def timestep(self,i):
        
        """
        Performs one time step within the orbit simulation. 
        The parameter [i] refers to which timestep we are currently at,
        meaning that 0 <= i < NSteps - 1.

        i       : int | Index of current step

        returns : void
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

        returns : void
        """

        # Loop
        for i in range(0,self.NSteps-1):
            # Performing the step
            self.timestep(i)     

        # Returning our arrays for plotting. We transpose the arrays again before returning,
        # which gives us the arrays with dimentions (2,NumPlanets,NSteps), letting us plot for all timesteps.
        return(self.r.T,self.v.T,self.a.T,self.t)
    
class NumericalOrbitFunction:

    def __init__(self,Filepath):

        """
        This class is a representation of the data stored from the numerical data.
        Objects of this class can be called to get the position of a given planet at a given time.

        Filepath : String | The path to the .npz file the class reads from.

        returns  : self
        """

        # Loading the config and r arrays from file
        npz = np.load(Filepath)

        # Setting total time, delta time, and number of time steps, from the read file
        self.config       = npz["config"]
        self.TotalTime    = self.config[1]
        self.dt           = self.config[2]
        self.NumSteps     = self.config[3]

        self.OrbitTimes   = npz["OrbitTimes"]

        # Setting r from read file
        self.r = npz["r"]

        # Colors we display the different planets with
        self.colors = [[0,0,1], [0.3,0,1], [0.4,0,1], [0.5,0,1], [0.6,0,1], [0.7,0,1], [0.8,0,1]]
        self.primary = [0.5,0,1]


    def __call__(self,t,p):

        """
        Method that returns the position of a given planet along the x and y axes at a given time.

        t       : float        | the desired point in time
        p       : int          | the desired planet index

        returns : Array(float) | the position of the given planet at the given time
        """

        if(t < 0):
            t = self.RotationTime - t

        # Finding the index of the given time
        Index = int(np.floor((t/self.TotalTime)*self.NumSteps))

        # Finding x and y positions at this index
        x = (self.r[0][p][Index])
        y = (self.r[1][p][Index])

        # Returning vector
        return(np.array([x,y]))
    
    def range(self,t0,t1):

        """
        Method that returns the positions of the planets along the x and y axes between
        two given points in time.

        t0       : float       | the desired starting time
        t1       : float       | the desired ending time

        returns : Array(float) | the positions of all planets between t0 and t1.
        """

        if(t0 < 0):
            t0 = self.RotationTime + t0
            t1 += self.RotationTime
        if(t1 < 0):
            t1 = self.RotationTime + t1
            t0 += self.RotationTime

        # Finding the indexes in the array for t0 and t1
        Index0 = int(np.floor((t0/self.TotalTime)*self.NumSteps))
        Index1 = int(np.floor((t1/self.TotalTime)*self.NumSteps))

        # Creating some empty lists
        x = []
        y = []

        # Filling information for x and y axes
        for p in range(len(self.r[0])):
            
            # Using list slicing to get all points at once [start:end]
            x.append(self.r[0][p][Index0:Index1])
            y.append(self.r[1][p][Index0:Index1])

        # Returning array
        return(np.array([x,y]))
    
    def GetColors(self):

        """
        Method that returns the colors of the planet orbits in the correct order.
        
        returns : list | list of planet colors in correct order
        """

        return([self.colors[0],self.colors[1],self.colors[4],self.colors[3],self.colors[5],self.colors[6],self.colors[2]])


if __name__ == "__main__":

    # Initializing AST
    Seed = utils.get_seed('bmthune')
    mission = SpaceMission(Seed)
    system = mission.system

    # Getting initial conditions
    R0 = system.initial_positions
    V0 = system.initial_velocities
    # Calculating the time the simulation will run. Here we assume that the orbit is a perfect circle, which it isn't, but it's very close.
    # To make sure we pass the 20 rotations mark, we multiply the time with 2
    OrbitTimes = (2*np.pi)*((system.semi_major_axes**3)/(const.G_sol*system.star_mass + system.masses))**(1/2)#np.linalg.norm(R0.T[0]) * 2 * np.pi/np.linalg.norm(V0.T[0])
    TotalTime = OrbitTimes[0] * 20 * 2

    # Instantiating the Numerical Orbit class (and running the loop)
    Orbit = NumericalOrbit(mission = mission,const = const, TotalTime = TotalTime, StepsPerYear = 10000, InitialPos = R0, InitialVel = V0)
    r,v,a,t = Orbit.loop()

    # Initializing plotting
    fig, ax = plt.subplots()

    # Plotting the orbits of the planets
    colors = Orbit.GetColors()
    for i in range(Orbit.NumPlanets):
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

    # Saving the array r, along with the total time, delta time, and number of timesteps.
    config = np.array([Orbit.T,Orbit.dt,Orbit.NSteps])
    np.savez("NumericalOrbitData",r = r,config = config,OrbitTimes = OrbitTimes)

    

    # TESTING ZONE!!

    eps = 1e-3
    relative_eps = 0.01
    TestCount = 0

    o_A = AnalyticalOrbit(SemiMajors = system.semi_major_axes,Eccentricities = system.eccentricities,AphelionAngles=system.aphelion_angles)
    r_As = o_A.Loop()

    for i in range(Orbit.NumPlanets):
        # First testing if final position is correct for all planets

        LastOrbitTime = (TotalTime % OrbitTimes[i])
        OrbitRatio = LastOrbitTime/OrbitTimes[i]
        Angle = 2*np.pi*OrbitRatio

        # Getting the angle of the final position
        theta = np.arctan(r[1][i][-1]/r[0][i][-1])
        if(r[0][i][-1] < 0):
            theta += np.pi

        # Getting the analytical radius for this angle
        r_N = np.linalg.norm(np.array([r[0][i][-1],r[1][i][-1]]))
        #r_A = r_As[int(np.floor(len(o_A.Angles)*OrbitRatio))]
        r_A = (system.semi_major_axes[i]*(1-system.eccentricities[i]**2))/(1+system.eccentricities[i]*np.cos(theta - (np.pi + system.aphelion_angles[i])))
        # Checking if this matches the numerical radius
        relative_error = abs(r_N - r_A)/r_A

        #print(r_N,r_A,Angle,np.linalg.norm(np.array([r[0][i][0],r[1][i][0]])))

        if not (relative_error < relative_eps):
            raise ValueError(f"Final radius for planet {i + 1} not correct!\nValue is {r_N} AU, but should be {r_A} AU, ratio is {relative_error}")
        TestCount += 1

        r_0   = np.linalg.norm(np.array([r[0][i][0],r[1][i][0]]))
        r_max = r_0
        r_min = r_0
        t_o   = 0

        # Collecting orbit time, max radius, and min radius
        for t in range(Orbit.NSteps):
            r_t = np.linalg.norm(np.array([r[0][i][t],r[1][i][t]]))

            if(r_t > r_max):
                r_max = r_t
            if(r_t < r_min):
                r_min = r_t

            if(t == 0):
                if(abs(r_0 - r_t) < eps):
                    t_o = Orbit.dt * t

        t_A = (2*np.pi*(system.semi_major_axes[i]**3/(const.G_sol*system.star_mass + system.masses[i])))

        if not (t_o - t_A < eps):
            raise ValueError(f"Final orbit time for planet {i + 1} not correct!\nValue is {t_o} AU, but should be {t_A} AU")
        

