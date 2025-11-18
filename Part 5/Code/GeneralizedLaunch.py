# Bruker ikke kodemal!!!
# Skrevet av Bastian Eggum Huuse og Bendik Thune

import matplotlib.pyplot as plt
import numpy as np
import pickle as pckl
from numba import njit
from PIL import Image

import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission


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
        self.TotalTime    = self.config[0]
        self.dt           = self.config[1]
        self.NumSteps     = int(self.config[2])

        self.OrbitTimes   = npz["OrbitTimes"]

        # Setting r from read file
        self.r = npz["r"]
        self.v = npz["v"]
        self.a = npz["a"]

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

        # Wrapping the t-value
        # if t is less than zero, we make it wrap around to the end of the simulation
        # This stops index-issues.
        if(t < 0):
            t = self.OrbitTimes[p] - t

        # Finding the index of the given time
        # This deserves an explanation. Since the positions are stored in an array with a length of NumSteps,
        # we can't just insert t into this array to get the value (since t is a floating number)
        # t/self.TotalTime gives us the percentage of the simulation the time t is at.
        # (if t/self.TotalTime = 0.5, t is halfway through the simulation).
        # We multiply this number with the total amount of steps, to get the closes time index to our current time.
        # We then floor that index (round down) and turn it into an integer.
        Index = int(np.floor(((t)/self.TotalTime)*self.NumSteps))

        # Finding x and y positions at this index
        x = (self.r[0][p][Index])
        y = (self.r[1][p][Index])

        # Returning vector
        return(np.array([x,y]))
    
    def GetVelocity(self,t,p):
        """
        Method that returns the velocity of a given planet along the x and y axes at a given time.

        t       : float        | the desired point in time
        p       : int          | the desired planet index

        returns : Array(float) | the velocity of the given planet at the given time
        """

        # Wrapping the t-value
        # if t is less than zero, we make it wrap around to the end of the simulation
        # This stops index-issues.
        if(t < 0):
            t = self.OrbitTimes[p] - t

        Index = int(np.floor((t/self.TotalTime)*self.NumSteps))

        # Finding x and y positions at this index
        x = (self.v[0][p][Index])
        y = (self.v[1][p][Index])

        # Returning vector
        return(np.array([x,y]))

    def range(self,t0,t1):

        """
        Method that returns the positions of the planets along the x and y axes between
        two given points in time.

        t0       : float        | the desired starting time
        t1       : float        | the desired ending time

        returns  : Array(float) | the positions of all planets between t0 and t1.
        """

        if(t0 < 0):
            t0 = self.OrbitTimes[0] + t0
            t1 += self.OrbitTimes[0]
        if(t1 < 0):
            t1 = self.OrbitTimes[0] + t1
            t0 += self.OrbitTimes[0]

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
    
class FuelRocket:

    def __init__(self,FuelMass,SpeedBoost,NumMotors,mission,NumParticles = 10**5,dt = 10**(-3) ):

        # Parameters
        self.mission = mission
        self.FuelMass = FuelMass
        self.RocketMass = self.mission.spacecraft_mass
        self.TotalMass = FuelMass + self.RocketMass
        self.SpeedBoost = SpeedBoost
        self.NumMotors = NumMotors
        self.dt = dt
        self.t = 0

        # We only care about velocity
        self.Velocity = 0
        self.counter = 0

        # Motor Parameters
        self.MotorLength = 10**(-6)
        self.NozzleLength = self.MotorLength * 0.25
        self.Temp = 3000
        self.NumParticles = NumParticles

        self.Thrust, self.FuelConsumption = self.SimulateEngine()

        # Lists used to save data for plotting later
        self.VelocityList = []
        self.FuelList = []

    def SimulateEngine(self):

        # Initializing motor
        print("Initializing motor...")

        # Calculating the force and fuelconsumption of one motor
        Force = 1.26686e-10
        FuelConsumption = 2.95044e-14

        # We have N motors, so to simulate this, we multiply by NumMotors
        # We assume that all the motors operate the same.
        TotalForce = Force * self.NumMotors
        TotalFuelConsumption = FuelConsumption * self.NumMotors

        #print(f"Calculated Force per motor : {Force:.5e}, Calculated Fuel Consumption per motor : {FuelConsumption:.5e}")
        #print(f"Calculated Force           : {TotalForce:.5e}, Calculated Fuel Consumption      : {TotalFuelConsumption:.5e}")
        return(TotalForce,TotalFuelConsumption)

    def TimeStep(self):

        # Eulering Velocity, Fuel, and time
        self.Velocity += ((self.Thrust/self.TotalMass)) * self.dt
        self.FuelMass -= self.FuelConsumption * self.dt
        self.TotalMass = self.FuelMass + self.RocketMass
        self.counter += 1
        self.t += self.dt

        # Adding to lists (Comment out to increase performance)
        self.VelocityList.append(self.Velocity)
        self.FuelList.append(self.FuelMass)

        # Printing info every 1000th timestep (uncomment to see info)
        if self.counter % 1000 == 0:
            pass
            #print(f"Current Velocity : {self.Velocity:.3f}, Current Fuel Mass : {self.FuelMass:.3f},Current time in seconds : {self.t:.1f}, Current time in minutes : {self.t/60:.1f}")

    def TimeLoop(self):

        while(self.Velocity < self.SpeedBoost):
            self.TimeStep()
            #break

            if(self.FuelMass <= 0):
                print("HOUSTON WE HAVE A PROBLEM.... \nBAAANG")
                break

class SimulationRocket(FuelRocket):

    def __init__(self,FuelMass,SpeedBoost,NumMotors,mission,NumParticles = 10**5,dt = 10**(-3),Graph = False):
        super().__init__(FuelMass,SpeedBoost,NumMotors,mission,NumParticles,dt) #Using class from FuelRocket 
        
        # Initializing variables
        self.Mission = mission
        self.System = self.Mission.system
        self.PlanetMass = self.System.masses[0] * const.m_sun
        self.PlanetRadius = self.System.radii[0] * 1000
        self.GravityConstant = const.G

        r_y = self.PlanetRadius
        self.v_x = ((2*np.pi)/(self.System.rotational_periods[0] * (86400))) * r_y
        # rotational_periods[0] is given in 24 hours, so we have to turn it into seconds.
        # 86400 is the amount of seconds in 24 hours
        
        # Position and Velocity are now vectors!!
        self.Position = np.array([0,r_y])
        self.Velocity = np.array([self.v_x,0.0])

        # Graphing lists
        self.Graph = Graph
        self.Positions = []
        self.Velocities = []

    def TimeStep(self):

        # Calculating the total acceleration of the vessel without direction
        GravityAcceleration = -self.GravityConstant*(self.PlanetMass/(np.linalg.norm(self.Position)**2))
        ThrustAcceleration = self.Thrust/self.TotalMass
        TotalAcceleration = ThrustAcceleration + GravityAcceleration

        # Adding a direction to the acceleration
        AccelerationDirection = self.Position/np.linalg.norm(self.Position)
        AccelerationVector = TotalAcceleration * AccelerationDirection

        # Tracking stats
        if self.Graph: 
            self.Positions.append(self.Position.copy())
            self.Velocities.append(self.Velocity.copy())

        # Eulering Velocity, Position, Fuel, and Time
        self.Velocity += AccelerationVector * self.dt
        self.Position += self.Velocity * self.dt
        self.FuelMass -= self.FuelConsumption * self.dt
        self.TotalMass = self.FuelMass + self.RocketMass
        self.t += self.dt

        # Note! Technically this rocket can exist INSIDE the planet, particularily at the beginning of the simulation
        # At the beginning the mass is (for some parameters) too high for the thrust to overtake the gravitational force,
        # causing the rocket to accelerate into the planet. This changes very little of the rest of the simulation
        # and also it's lowkey kinda funny, so we have decided not to write code to change this. But if we would change we
        # would just set a hard boundary.

        # Calculating new Escape velocity
        self.SpeedBoost = ((2 *(self.PlanetMass)*self.GravityConstant) / (np.linalg.norm(self.Position)))**(1/2)

    def TimeLoop(self):
        while(np.linalg.norm(self.Velocity + np.array([self.v_x,0])) < self.SpeedBoost):
            
            self.TimeStep()

            if(self.FuelMass <= 0):
                print("HOUSTON WE HAVE A PROBLEM.... \nBAAANG")
                print("     _.-^^---....,,-- \n _--                  --_\n<                        >)\n|                         |\n \\._                   _./\n    ```--. . , ; .--'''\n          | |   |\n       .-=||  | |=-.\n       `-=#$%&%$#=-'\n          | ;  :|\n _____.,-#%&$@%#&#~,._____ ")
                raise RuntimeError
                

    def StarPosition(self,Pos,Vel):
        
        # Our code already takes into account of the rotation of our planet
        # so we don't have to take that into account when calculating new velocity

        # According to the image in the problem description, we launch along the x axis in the solar system frame
        # This means that our axes have to be switched!!

        Pos[0],Pos[1] = Pos[1],Pos[0]
        Vel[0],Vel[1] = Vel[1],Vel[0]

        # Adjusting the position and velocity to be in Astronomical units
        Pos *= (1/const.AU)
        Vel *= ((60*60*24*365)/const.AU)
        
        # Getting planet position (in AU)
        r_p = self.System.initial_positions[:,0]

        # Setting our position in the solar system frame 
        r = r_p + Pos

        # Getting planet velocity
        v_p = self.System.initial_velocities[:,0]

        # Uppdating Possision with plantets orbit speed
        r += v_p * self.t/(60*60*24*365)

        # Setting our velocity in the solar system frame
        v = v_p + Vel
        
        return r,v
    
class GeneralizedRocket(SimulationRocket):

    def __init__(self,mission,FuelMass,SpeedBoost,NumMotors,NumParticles = 10**5,dt = 10**(-3)):

        """
        Class that simulates a rocket launch from our planet, at a given time t and from a given angle theta

        Most functionality from this class is inherited from the class SimulationRocket, which was written for Part 1.
        Comments have only been added for new functionality
        """

        super().__init__(FuelMass,SpeedBoost,NumMotors, mission,NumParticles,dt)

        # Tracking some stats about the planet
        self.mission = mission
        self.system = self.mission.system

        self.FileName = "NumericalOrbitData.npz"
        self.R_planets = NumericalOrbitFunction(self.FileName)

        # Getting Orbit- and Rotationtime (both in years)
        self.OrbitTime = (2*np.pi)*((self.system.semi_major_axes[0]**3)/(const.G_sol*(self.system.star_mass + self.system.masses[0])))**(1/2)
        self.RotationTime = self.system.rotational_periods[0] / 365

    def WrapTime(self, t):

        """
        Method that wraps time.

        t       : float | time before wrapping

        returns : float | time after wrapping
        """

        # while time is below 0, add the time of 1 rotation
        while t < 0:
            t += self.OrbitTime

        # while time is above (or equal to) OrbitTime, subtract 1 rotation
        while t >= self.OrbitTime:
            t -= self.OrbitTime

        # This gives us 0 <= t < self.Orbittime
        return t

    def SolarSystemPosition(self,t_0,theta):

        """
        Method that calculates position after launch when launching at time t_0 and from angle theta

        t_0   : float | Time at beginning of launch in Years
        theta : float | Angle of launch site (around the equator) in radians

        returns :
        Array(float) | Position after launch (in solar system coordinates) in AU
        Array(float) | Velocity after launch (in solar system coordinates) in AU/Y
        """

        # Wrapping time
        t_0 = self.WrapTime(t_0)

        # Copying arrays so we don't mess with them
        Pos = self.Position.copy()
        Vel = self.Velocity.copy()

        # We have a position on the y axis, radially outward from the planet, and a position on the x-axis, normally along this axis.
        # Both of these are in meters, so we must change to AU:
        Pos *= (1/const.AU)
        Vel *= ((60*60*24*365)/const.AU)

        # Now we can calculate the velocity.
        # First we get the current velocity of the planet:
        v_p = self.R_planets.GetVelocity(t_0,0)
        # Then, we swap the velocity axes around
        Vel[0], Vel[1] = Vel[1], Vel[0]
        # And finally add the velocities together:
        SolarSystemVel = v_p + np.array([Vel[0] * np.cos(theta) - Vel[1] * np.sin(theta), Vel[0] * np.sin(theta) + Vel[1] * np.cos(theta)])

        # Now, we want to get the current planet position.
        r_p = self.R_planets(t_0,0)
        # Note that we swap the x and y axes here (just like the original program)
        Pos[0], Pos[1] = Pos[1], Pos[0]
        # We rotate the local positional vector by the angle theta (Notice that this is just multiplying with the rotation matrix), and add the solar system position of the planet.
        SolarSystemPos = r_p + np.array([Pos[0] * np.cos(theta) - Pos[1] * np.sin(theta), Pos[0] * np.sin(theta) + Pos[1] * np.cos(theta)])
        # Adding the contribution of the planet's orbit speed
        SolarSystemPos += v_p * (self.t / (60 * 60 * 24 * 365))

        r_hat = r_p / np.linalg.norm(r_p)
        # Finally, we also get the launch position for plotting later
        self.LaunchPos = r_p + ((self.system.radii[0] * 1000)/const.AU) * r_hat#np.array([np.cos(theta),np.sin(theta)])#r_hat

        # Now we can return these values:
        return SolarSystemPos,SolarSystemVel    
    
@njit
def Doppler2V(dLambda, Lambda_0):
    """
    Finds the radiall speed of one star by messuring the dopple effect
    Inputs:
    dLambda: Messured Doppler shift (nm)
    Lambda_0: Original spectral line(nm)
    Outputt: Velocity in radial direction gotten from dopplershift (m/s)
    """
    v_r = dLambda/Lambda_0 * const.c #the formula for dopler shift
    return v_r


@njit
def RelVel(dLambda_1, dLambda_2,Lambda_0):
    """
    Creates a tupple with the velocity relativ to two stars
    Inputs:
    dLambda_1: Messured dopplereffect of star 1 (nm)
    dLambda_2: Messured dopplereffect of star 2 (nm)
    Lambda_0: Original spectral line (nm)
    Outputt: Vellocity given in the basis of star 1 and 2  (m/s)
    """

    v_1 = Doppler2V(dLambda_1,Lambda_0) 
    v_2 = Doppler2V(dLambda_2, Lambda_0)
    return (v_1,v_2)

@njit
def BasisShift(phi_1,phi_2, v):
    """
    Changes basis from star basis to xy
    Inputs:
    phi_1: angle to star 1
    phi_2: angle to star 2
    v: velocity given in basis of speed away from star 1 and 2
    Outputt: velocity in xy basis
    """
    #Here we change the basis from stars to x,y notice adding pi, this is because the angles gives us the direction toward the stars
    #while the basis goes the other way
    v_x = np.cos(phi_1 + np.pi)*v[0] + np.cos(phi_2+np.pi) *v[1]
    v_y = np.sin(phi_1 + np.pi)*v[0] + np.sin(phi_2 +np.pi) *v[1]
    return v_x, v_y

@njit
def RockVel(v_r,v_s):
    """
    Find the relativ velocity of rocket to sun
    Inputs:
    v_r: vellocity of rocket compared to stars 
    v_s: vellocity of sun compared to stars (v_s and v_r must have same unit)
    outputt: Vellocity of rocket compared to sun (same unit as sent in)
    """
    v_x = v_r[0] - v_s[0]
    v_y = v_r[1] - v_s[1]
    return v_x, v_y #not nesesarly in xy basis

@njit
def DopplerVelocity(dLambda, Lambda_0, phi_1, phi_2):
    """
    Runns all functions above to find the vellocity 
    Inputs:
    dLambda: an array of length 4 that has the dopplershifts from each star on the sun and rocket (in nm)
    Lambda_0: Original spectral line in nm
    phi_1: angle to star 1
    phi_2: angle to star 2
    Outputt: Velocity of rocket in xy basis compared to sun (AU/Y)
    """
    #first i find the velocity of rocket at sun in star basis
    v_s = RelVel(dLambda[0],dLambda[1],Lambda_0)
    v_r = RelVel(dLambda[2],dLambda[3],Lambda_0)
    
    #Find the velocity of rocket compared to sun in star basis
    V = RockVel(v_r, v_s)
    
    #change the basis to xy basis
    V_r = BasisShift(phi_1,phi_2, V)
    
    #returns the velocity converted to AU/Y
    return V_r[0]*60*60*24*365/const.AU, V_r[1]*60*60*24*365/const.AU

@njit
def CompareImages(Image1, Image2):

    '''
    Method that compares two images by taking the sum of the squares of the differences of their RGB values

    Parameters:
    Image1 : Array(float) | the first image of the two we want to compare, converted to an array
    Image2 : Array(float) | the second image of the two we want to compare, converted to an array

    returns : float       | sum of the squares of the differences of the RGB values of the two images
    '''

    # Copying arrays so we don't modify the originals
    Img1 = Image1.copy().astype("float32")
    Img2 = Image2.copy().astype("float32")

    Diff = 0

    # Comparing iteratively (Going through each pixel)
    for i in range(len(Img1)):
        for j in range(len(Img1[i])):

            # Computing the differences of the R, G, and B values
            dR = Img1[i][j][0] - Img2[i][j][0]
            dG = Img1[i][j][1] - Img2[i][j][1]
            dB = Img1[i][j][2] - Img2[i][j][2]

            # Summing the squares of these
            dC = (dR**2 + dG**2 + dB**2)
            
            # Adding this sum to the total difference
            Diff += dC

    return Diff

@njit
def CompareImageRange(InputImage,SkyImageData,N = 360):

    '''
    Comparing an input image with 360 other images, one for each degree.
    The method then finds the image with the lowest total difference, which is the closest match

    Parameters:
    InputImage   : Array(float) | the image we want to approximate, converted to an array
    SkyImageData : Array(float) | an array containing 360 images, one for each degree (0-359)
    N            : int          | how many of the degrees we want to compare InputImage to (0 <= N < 360)
    '''

    # Getting the N first images
    ComparisonImages = SkyImageData[:N]
    Diffs = np.zeros(N)

    # Calculating all differences
    for i in range(len(Diffs)):
        Diffs[i] = CompareImages(InputImage,ComparisonImages[i])

    # Finding the lowest difference (index)
    MinIndex = 0
    for i in range(1,len(Diffs)):

        if(Diffs[i] < Diffs[MinIndex]):
            MinIndex = i

    # Since we have 360 different angles indexed,
    # we can return the angle by simply returning the index
    # (Note that this gives angle in degrees)

    return MinIndex

@njit
def TupleDiff(Tuple1,Tuple2):

    """
    Support function that returns the difference between two tuples containing float values
    (Note that both tuple functions also work with arrays)

    Parameters:
    Tuple1 : tuple | The first tuple
    Tuple2 : tuple | The second duple

    Returns:
    tuple | The difference between the two tuples
    """

    Coord0 = Tuple1[0] - Tuple2[0]
    Coord1 = Tuple1[1] - Tuple2[1]
    return Coord0,Coord1

@njit
def TupleNorm(Tuple):

    """
    Support function that returns the norm of a tuple containing two float values

    Parameters:
    Tuple : tuple | The tuple

    Returns:
    float | The norm of the tuple
    """

    return (Tuple[0]**2 + Tuple[1]**2)**(1/2)

@njit
def PlanetDifference(Pos,PlanetPositions,TruePlanetDistances,printer = True):

    """
    Method that computes the difference between the distances from the planets to a given point, and what those distances should be.
    
    Parameters:
    Pos                 : Array(float) | The current position we want to check
    PlanetPositions     : Array(float) | List of planet positions
    TruePlanetDistances : Array(float) | List of expected distances

    returns: float | The difference between computed distances and expected distances (the sum of these)
    """

    PlanetDeviations = np.zeros((len(PlanetPositions),2))

    # Finding the deviation vector between the given point and all planets
    for i in range(len(PlanetDeviations)):
        Deviation = TupleDiff(Pos,PlanetPositions[i])
        PlanetDeviations[i][0] = Deviation[0]
        PlanetDeviations[i][1] = Deviation[1]

    TotalDifference = 0

    # Computing the sum of the difference between each computed distance, and each expected distance
    for i in range(len(PlanetDeviations)):

        TotalDifference += abs(TupleNorm(PlanetDeviations[i]) - TruePlanetDistances[i])

    if(printer):
        print("Total: ",TotalDifference)

    return TotalDifference

@njit
def BinaryLeastSquares(StarDistance,Range,PlanetPositions,TruePlanetDistances):

    """
    Algorithm that finds the point where expected distances match computed distances (closest). This would then be the point in space where the rocket is.
    This point is found through a combination binary search and least squares approach.

    We create a circle with center in the star, and radius equal to the distance the rocket is from the sun.
    We place a point on each side of this circle, and find the point closest to the expected position through least squares
    We then select the half of the circle the closest point is on, and then call the method recursively on that part of the circle.

    Parameters:
    StarDistance        : float        | The distance between the rocket and the star
    Range               : float        | The range of angles we are currently checking (Begins as 0-2pi)
    PlanetPositions     : Array(float) | The positions of all the planets
    TruePlanetDistances : Array(float) | The expected distances of each planet 
    
    returns: float | The angle from the star the position of the rocket is located (in radians)
    """

    # If the range of angles has a length of 1, we simply return that angle (We've found the closest angle)
    if len(Range) == 1:
        return Range[0]

    # Splitting Range in two
    HalfLength = int(len(Range)/2)
    Range1 = np.zeros(HalfLength)
    Range2 = np.zeros(HalfLength)

    # Looping through and duplicating the half ranges
    for i in np.arange(HalfLength*2):
        
        if i < HalfLength:
            Range1[i] = Range[i]
        else:
            Range2[i - HalfLength] = Range[i]

    # Finding the Positions we want to check
    HalfIndex = int(HalfLength/2)
    Position1 = (StarDistance * np.cos(Range1[HalfIndex]),StarDistance * np.sin(Range1[HalfIndex]))
    Position2 = (StarDistance * np.cos(Range2[HalfIndex]),StarDistance * np.sin(Range2[HalfIndex]))

    # Finding the least squares difference for both of these points
    Diff1 = PlanetDifference(Position1,PlanetPositions,TruePlanetDistances)
    Diff2 = PlanetDifference(Position2,PlanetPositions,TruePlanetDistances)

    # Recursing with the smallest one
    if Diff1 <= Diff2:
        Angle = BinaryLeastSquares(StarDistance,Range1,PlanetPositions,TruePlanetDistances)
    else:
        Angle = BinaryLeastSquares(StarDistance,Range2,PlanetPositions,TruePlanetDistances)

    return Angle

def TrilaterationAlgorithm(t,Distances):
    
    """
    Main function that actually calls the algorithm.

    Parameters:
    t         : float        | Time of launch (in years)
    Distances : Array(float) | Array containing the Expected distances of all planets (and the star)

    returns: Array(float) | The position of the rocket after launch (in AU)
    """

    # Instantiating our Position function
    FilePath = "NumericalOrbitData.npz"
    PlanetPositionFunction = NumericalOrbitFunction(FilePath)

    # Defining star-distance and the range of angles
    StarDistance = Distances[-1]
    Grain = 15
    Range = np.linspace(0,2*np.pi,2**(Grain))
    
    # Defining an array of distances that doesn't include the sun
    PlanetDistances = np.zeros((len(Distances)-1))
    for i in range(len(PlanetDistances)):
        PlanetDistances[i] = Distances[i]

    # Defining an array of positions, one for each planet.
    PlanetPositions = np.zeros((len(PlanetDistances),2))
    for p in range(len(PlanetPositions)):
        PlanetPositions[p] = PlanetPositionFunction(t,p)

    # Running the algorithm to find the angle from the sun at which our position is at.
    Angle = BinaryLeastSquares(StarDistance,Range,PlanetPositions,PlanetDistances)

    # We want to return the position, instead of just the angle
    return (StarDistance * np.cos(Angle),StarDistance * np.sin(Angle))

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
        self.primary = [0.0,0,1]

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

def main(mission, t_0):

    
    """
    This function does three things.

    First, it creates the file containing the planet positions
    """
    # system = mission.system
    # # Getting initial conditions
    # R0 = system.initial_positions
    # V0 = system.initial_velocities
    # # Calculating the time the simulation will run. Here we assume that the orbit is a perfect circle, which it isn't, but it's very close.
    # # To make sure we pass the 20 rotations mark, we multiply the time with 2
    # OrbitTimes = (2*np.pi)*((system.semi_major_axes**3)/(const.G_sol*(system.star_mass + system.masses)))**(1/2)#np.linalg.norm(R0.T[0]) * 2 * np.pi/np.linalg.norm(V0.T[0])
    # TotalTime = OrbitTimes[0] * 20 * 2

    # # Instantiating the Numerical Orbit class (and running the loop)
    # # We found that 10000 steps per year is sufficient, as all tests provide reasonable results with these parameters
    # # Increasing the steps per year would then only reduce performance.
    # Orbit = NumericalOrbit(mission = mission,const = const, TotalTime = TotalTime, StepsPerYear = 10000, InitialPos = R0, InitialVel = V0)
    # r,v,a,t = Orbit.loop()

    # # Saving the array r, the planets' orbit times, and the total time, delta time, and number of timesteps.
    # # We use this file later, to circumvent having to run the simulation again.
    # config = np.array([Orbit.T,Orbit.dt,Orbit.NSteps])
    # np.savez("NumericalOrbitData",r = r, v = v, a = a,config = config,OrbitTimes = OrbitTimes)

    """
    Then, it simulates the launch of the rocket
    """

    # Creating rocket instance
    NumMotors = int((1000000**3)/60) # 1/10 qube meter grid :)
    Fuel = 190000
    Particles = 10**5
    EscapeVelocity = ((2 *(mission.system.masses[0]*const.m_sun)*const.G) / (mission.system.radii[0] * 1000))**(1/2)

    # Creating a generalized rocket
    GenRocket = GeneralizedRocket(mission=mission,FuelMass=Fuel,SpeedBoost=EscapeVelocity,NumMotors=NumMotors,NumParticles=Particles)

    # Looping rocket
    GenRocket.TimeLoop()

    print(f"Fuel after launch: {GenRocket.FuelMass}")

    t_0 = t_0
    Sim_Pos0, Sim_Vel0 = GenRocket.Position.copy(),GenRocket.Velocity.copy()
    Sim_Pos, Sim_Vel = GenRocket.StarPosition(Sim_Pos0,Sim_Vel0)

    r_p = GenRocket.R_planets(t_0,0)
    z = r_p[0] + 1j*r_p[1]
    Gen_Pos, Gen_Vel = GenRocket.SolarSystemPosition(t_0, np.angle(z))

    print(f"\nGeneralized Position in solar system frame at t = {t_0}: [x : {Gen_Pos[0]:.3f} AU, y : {Gen_Pos[1]:.2e} AU]")
    print(f"Generalized Velocity in solar system frame at t = {t_0}: [x : {Gen_Vel[0]:.3f} AU/Y, y : {Gen_Vel[1]:.3f} AU/Y]")

    print(f"\nSpecialized Position in solar system frame at t = {t_0}: [x : {Sim_Pos[0]:.3f} AU, y : {Sim_Pos[1]:.2e} AU]")
    print(f"Specialized Velocity in solar system frame at t = {t_0}: [x : {Sim_Vel[0]:.3f} AU/Y, y : {Sim_Vel[1]:.3f} AU/Y]")

    print()

    mission.set_launch_parameters(
        thrust = GenRocket.Thrust,
        mass_loss_rate = GenRocket.FuelConsumption,
        initial_fuel_mass = Fuel,
        estimated_launch_duration = GenRocket.t + 1,
        launch_position = GenRocket.LaunchPos,
        time_of_launch =   t_0 -1*GenRocket.R_planets.dt
        )
    
    mission.launch_rocket(10**(-3))
    
    mission.verify_launch_result(Gen_Pos)

    """
    Then, it performs a manual orientation.
    """

    t_1 = t_0 + (GenRocket.t/(60*60*24*365))

    # Collecting information:
    # Rotation:
    InputPath = "Sky_Image.png"
    mission.take_picture(InputPath,"himmelkule.npy")
    InputImage = Image.open(InputPath)
    InputImageArray = np.array(InputImage)
    SkyImageData = np.load("SkyImageData.npy")
    # Velocity:
    lambda_0 = 656.3 #Hydrogen spectral line
    lambda_1, lambda_2 = mission.star_doppler_shifts_at_sun
    lambda_3, lambda_4 = mission.measure_star_doppler_shifts()
    dlambda = (lambda_1,lambda_2,lambda_3,lambda_4)
    phi_1, phi_2 = mission.star_direction_angles
    phi_1, phi_2 = (np.deg2rad(phi_1), np.deg2rad(phi_2))
    # Position
    Distances = mission.measure_distances()

    # Performing computations
    Orientation_Rotation = CompareImageRange(InputImageArray,SkyImageData,N = 360)
    Orientation_Velocity = DopplerVelocity(dlambda, lambda_0, phi_1, phi_2)
    Orientation_Position = TrilaterationAlgorithm(t_1,Distances)
    print(f"Distances {Distances}")
    print(Orientation_Position)

    # Verifying results
    mission.verify_manual_orientation(Orientation_Position,Orientation_Velocity,Orientation_Rotation)

if __name__ == "__main__":
    seed = utils.get_seed('bmthune')
    mission = SpaceMission(seed)  

    main(mission, t_0= 0.2*11)

    with open ("Mission.pkl", 'wb') as file:
        pckl.dump(mission, file)


"""
String that runs code: python GeneralizedLaunch.py

Initializing motor...
Calculated Force per motor : 1.26686e-10, Calculated Fuel Consumption per motor : 2.95044e-14
Calculated Force           : 2.11143e+06, Calculated Fuel Consumption      : 4.91740e+02

Generalized Position in solar system frame at t = 2.8: [x : 2.577 AU, y : -1.12e+00 AU]
Generalized Velocity in solar system frame at t = 2.8: [x : 4.647 AU/Y, y : 4.373 AU/Y]

Specialized Position in solar system frame at t = 2.8: [x : 2.814 AU, y : 6.98e-05 AU]
Specialized Velocity in solar system frame at t = 2.8: [x : 2.480 AU/Y, y : 5.856 AU/Y]

Rocket was moved down by 423.957 m to stand on planet surface.
New launch parameters set.
Launch completed, reached escape velocity in 379.294 s.
Your spacecraft position was satisfyingly calculated. Well done!
*** Achievement unlocked: No free launch! ***
Picture written to Sky_Image.png.
Pointing angle after launch correctly calculated. Well done!
Velocity after launch correctly calculated. Well done!
Position after launch correctly calculated. Well done!
Your manually inferred orientation was satisfyingly calculated. Well done!
*** Achievement unlocked: Well-oriented! ***
"""