# BRUKER IKKE KODEMAL!!!!
# Skrevet av Bastian og Bendik (BaBe space center)

# Imports
import  numpy        as     np
import  matplotlib.pyplot as plt
# AST imports
import ast2000tools.constants as const
import ast2000tools.utils as utils
# From imports
from .FuelRocket import FuelRocket
from ast2000tools.space_mission import SpaceMission

# Ast init
seed = utils.get_seed('bmthune')
mission = SpaceMission(seed)


class SimulationRocket(FuelRocket):

    def __init__(self,FuelMass,SpeedBoost,NumMotors,NumParticles = 10**5,dt = 10**(-3),Graph = False):
        super().__init__(FuelMass,SpeedBoost,NumMotors,NumParticles,dt) #Using class from FuelRocket 
        
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

if __name__ == "__main__":
    

    # Creating rocket instance
    NumMotors = int((1000000**3)/60) # 1/10 qube meter grid :)
    Fuel = 190000
    Particles = 10**5
    EscapeVelocity = ((2 *(mission.system.masses[0]*const.m_sun)*const.G) / (mission.system.radii[0] * 1000))**(1/2)
    TrackValues = True

    TestRocket = SimulationRocket(FuelMass=Fuel,SpeedBoost=EscapeVelocity,NumMotors=NumMotors,NumParticles=Particles,Graph = TrackValues)
    TestRocket.TimeLoop()

    print(f"\nThe rocket has reached the escape velocity of {TestRocket.SpeedBoost:.2f} m/s!!! (or crashed)\n")
    
    # print(f"Info time :)")
    # print(f"Position post launch : {f"[x : {TestRocket.Position[0]:.2f} m, y : {TestRocket.Position[1]:.2f} m]":>40}")
    # print(f"Velocity post launch : {f"[x : {TestRocket.Velocity[0]:.2f} m/s, y : {TestRocket.Velocity[1]:.2f} m/s]":>40}")
    # print(f"Total time of launch : {f"[Seconds : {TestRocket.t:.1f} s, Minutes : {TestRocket.t/60:.1f} min]":>40}")
    # print(f"Total Fuel Consumed  : {f"{Fuel - TestRocket.FuelMass:.2f} kg":>40}")
    # print(f"Total Fuel Remaining : {f"{TestRocket.FuelMass:.2f} kg":>40}")
    # print(f"Total Rocket Mass    : {f"{TestRocket.TotalMass:.2f} kg":>40}\n")

    # Plotting
    fig, ax = plt.subplots() # Plotting init

    # Adding our planet
    planet = plt.Circle((0, 0), mission.system.radii[0] * 1000, color='r')
    ax.add_patch(planet)

    # Plotting our graph
    Points = list(zip(*TestRocket.Positions))
    N_Points = len(Points[0])
    ax.plot(Points[0][:N_Points],Points[1][:N_Points])
    plt.axis('equal')
    plt.show()

    r,v = TestRocket.StarPosition(TestRocket.Position,TestRocket.Velocity)
    print(f"Position in solar system frame : [x_f : {TestRocket.Position[0]}, x : {r[0]:.3f} AU, y : {r[1]:.2e} AU]")
    print(f"Velocity in solar system frame : [x : {v[0]:.3f} AU/Y, y : {v[1]:.3f} AU/Y]")

    mission.set_launch_parameters(
        thrust = TestRocket.Thrust,
        mass_loss_rate = TestRocket.FuelConsumption,
        initial_fuel_mass = Fuel,
        estimated_launch_duration = TestRocket.t + 1,
        launch_position = mission.system.initial_positions[:,0] + np.array([(mission.system.radii[0]*1000)/const.AU,0]),
        time_of_launch = 0
        )
    
    mission.launch_rocket(10**(-3))
    
    mission.verify_launch_result(r)
    