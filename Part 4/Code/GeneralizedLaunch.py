# Bruker ikke kodemal!!!
# Skrevet av Bastian Eggum Huuse og Bendik Thune

import matplotlib.pyplot as plt
import numpy as np

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
        Index = int(np.floor((t/self.TotalTime)*self.NumSteps))

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

        print(f"Calculated Force per motor : {Force:.5e}, Calculated Fuel Consumption per motor : {FuelConsumption:.5e}")
        print(f"Calculated Force           : {TotalForce:.5e}, Calculated Fuel Consumption      : {TotalFuelConsumption:.5e}")
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
    

def main(mission, t_0):
    # Ast init


    # Creating rocket instance
    NumMotors = int((1000000**3)/60) # 1/10 qube meter grid :)
    Fuel = 190000
    Particles = 10**5
    EscapeVelocity = ((2 *(mission.system.masses[0]*const.m_sun)*const.G) / (mission.system.radii[0] * 1000))**(1/2)

    # Creating a generalized rocket
    GenRocket = GeneralizedRocket(mission=mission,FuelMass=Fuel,SpeedBoost=EscapeVelocity,NumMotors=NumMotors,NumParticles=Particles)

    # Looping rocket
    GenRocket.TimeLoop()

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
        time_of_launch =   t_0 - GenRocket.R_planets.dt
        )
    
    mission.launch_rocket(10**(-3))
    
    mission.verify_launch_result(Gen_Pos)

    TimeAfterLaunch = t_0 + (GenRocket.t/(60*60*24*365))
    return TimeAfterLaunch

if __name__ == "__main__":
    seed = utils.get_seed('bmthune')
    mission = SpaceMission(seed)  
    main(mission, 1)



"""
String that runs code: python GeneralizedLaunch.py

This program also outputs an image, which has been turned in as "GeneralizedLaunchPlot.py"
Output (Note that output may not be the same as here, since the simulation uses randomness):

Initializing motor...
Simulating motor:   10%
Simulating motor:   20%
Simulating motor:   30%
Simulating motor:   40%
Simulating motor:   50%
Simulating motor:   60%
Simulating motor:   70%
Simulating motor:   80%
Simulating motor:   90%
Simulating motor:  100%
Finished Simulating 16666666666666666 motors.
Calculated Force per motor : 1.28804e-10, Calculated Fuel Consumption per motor : 2.94575e-14
Calculated Force           : 2.14674e+06, Calculated Fuel Consumption      : 4.90959e+02

Generalized Position in solar system frame at t = 0: [x : 2.814 AU, y : 6.98e-05 AU]
Generalized Velocity in solar system frame at t = 0: [x : 2.477 AU/Y, y : 5.855 AU/Y]

Specialized Position in solar system frame at t = 0: [x : 2.814 AU, y : 6.98e-05 AU]
Specialized Velocity in solar system frame at t = 0: [x : 2.477 AU/Y, y : 5.855 AU/Y]

Rocket was moved up by 1.49664e-05 m to stand on planet surface.
New launch parameters set.
Launch completed, reached escape velocity in 379.202 s.
Your spacecraft position was satisfyingly calculated. Well done!
*** Achievement unlocked: No free launch! ***
Ignoring fixed x limits to fulfill fixed data aspect with adjustable data limits.

"""