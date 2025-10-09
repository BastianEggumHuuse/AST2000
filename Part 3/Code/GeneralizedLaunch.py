import sys
sys.path.insert(0, "../..")

import matplotlib.pyplot as plt
import numpy as np

import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

from Del_1_Code.SimulationRocket import SimulationRocket
from OrbitPlotNumerical import NumericalOrbitFunction

class GeneralizedRocket(SimulationRocket):

    def __init__(self,mission,FuelMass,SpeedBoost,NumMotors,NumParticles = 10**5,dt = 10**(-3)):
        super().__init__(FuelMass,SpeedBoost,NumMotors,NumParticles,dt)

        # Tracking some stats about the planet
        self.mission = mission
        self.system = self.mission.system

        self.FileName = "NumericalOrbitData.npz"
        self.R_planets = NumericalOrbitFunction(self.FileName)

        # Getting Orbit- and Rotationtime (both in years)
        self.OrbitTime = (2*np.pi)*((self.system.semi_major_axes[0]**3)/(const.G_sol*(self.system.star_mass + self.system.masses[0])))**(1/2)
        self.RotationTime = self.system.rotational_periods[0] / 365
        

    def WrapTime(self, t):

        if t < 0:
            t = self.OrbitTime - t

        while t > self.OrbitTime:
            t -= self.OrbitTime

        return t

    def SolarSystemPosition(self,t_0,theta):

        t_0 = self.WrapTime(t_0)

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

        self.LaunchPos = r_p + ((self.system.radii[0] * 1000)/const.AU) * np.array([np.cos(theta),np.sin(theta)])

        # Now we can return these values:
        return SolarSystemPos,SolarSystemVel

    
if __name__ == "__main__":

    # Ast init
    seed = utils.get_seed('bmthune')
    mission = SpaceMission(seed)

    # Creating rocket instance
    NumMotors = int((1000000**3)/60) # 1/10 qube meter grid :)
    Fuel = 190000
    Particles = 10**5
    EscapeVelocity = ((2 *(mission.system.masses[0]*const.m_sun)*const.G) / (mission.system.radii[0] * 1000))**(1/2)

    # Creating a generalized rocket
    GenRocket = GeneralizedRocket(mission=mission,FuelMass=Fuel,SpeedBoost=EscapeVelocity,NumMotors=NumMotors,NumParticles=Particles)

    # Looping rocket
    GenRocket.TimeLoop()

    Sim_Pos0, Sim_Vel0 = GenRocket.Position.copy(),GenRocket.Velocity.copy()
    Sim_Pos, Sim_Vel = GenRocket.StarPosition(Sim_Pos0,Sim_Vel0)
    Gen_Pos, Gen_Vel = GenRocket.SolarSystemPosition(0,0)

    print(f"\nGeneralized Position in solar system frame at t = 0: [x : {Gen_Pos[0]:.3f} AU, y : {Gen_Pos[1]:.2e} AU]")
    print(f"Generalized Velocity in solar system frame at t = 0: [x : {Gen_Vel[0]:.3f} AU/Y, y : {Gen_Vel[1]:.3f} AU/Y]")

    print(f"\nSpecialized Position in solar system frame at t = 0: [x : {Sim_Pos[0]:.3f} AU, y : {Sim_Pos[1]:.2e} AU]")
    print(f"Specialized Velocity in solar system frame at t = 0: [x : {Sim_Vel[0]:.3f} AU/Y, y : {Sim_Vel[1]:.3f} AU/Y]")

    print()
    mission.set_launch_parameters(
        thrust = GenRocket.Thrust,
        mass_loss_rate = GenRocket.FuelConsumption,
        initial_fuel_mass = Fuel,
        estimated_launch_duration = GenRocket.t + 1,
        launch_position = mission.system.initial_positions[:,0] + np.array([(mission.system.radii[0]*1000)/const.AU,0]),
        time_of_launch = 0
        )
    
    mission.launch_rocket(10**(-3))
    
    mission.verify_launch_result(Gen_Pos)


    # Plotting
    t_0 = 1
    theta = np.pi/2

    # Getting info about Rocket and about planet
    Gen_Pos, Gen_Vel = GenRocket.SolarSystemPosition(t_0,theta)
    r_N = GenRocket.R_planets.range(0,t_0 + GenRocket.R_planets.dt)
    
    Planet_Pos_Pre = GenRocket.R_planets(t_0,0)
    Planet_Pos_Post = Planet_Pos_Pre + GenRocket.R_planets.GetVelocity(t_0,0) * (GenRocket.t / (60 * 60 * 24 * 365))
    Planet_Pos_Mid = (Planet_Pos_Pre + Planet_Pos_Post)/2
    Planet_r = (mission.system.radii[0] * 1000)/const.AU
    Radii = 2.5
    xlim = [Planet_Pos_Mid[0] - Planet_r * Radii, Planet_Pos_Mid[0] + Planet_r * Radii]
    ylim = [Planet_Pos_Mid[1] - Planet_r * Radii, Planet_Pos_Mid[1] + Planet_r * Radii]

    # Initializing plotting
    fig, ax = plt.subplots()

    # Plotting the orbit of the planet
    color = GenRocket.R_planets.primary
    ax.plot(r_N[0][0],r_N[1][0], color = color)

    # Plotting launch position and final position
    ax.plot(Gen_Pos[0],Gen_Pos[1],".",color = "Red",label = "Etter oppskytning")
    ax.plot(GenRocket.LaunchPos[0],GenRocket.LaunchPos[1],".",color = "Blue",label = "Før oppskytning")

    # Adding patches
    Planet_Pre = plt.Circle(Planet_Pos_Pre, Planet_r, color = 'Green',label = "Planet før oppskytning")
    Planet_Post = plt.Circle(Planet_Pos_Post, Planet_r, color = 'Lime',label = "Planet etter oppskytning")
    ax.add_patch(Planet_Pre)
    ax.add_patch(Planet_Post)

    # Defining limits
    plt.axis('equal')
    plt.xlim(xlim)
    plt.ylim(ylim)

    # Adding title, and axis labels
    plt.title(f"Oppskytning langs y-aksen ved t_0 = {t_0}")
    plt.xlabel("Posisjon langs x-aksen [AU]")
    plt.ylabel("Posisjon langs y-aksen [AU]")
    plt.legend(loc = "lower right")

    # showing plot
    plt.show()

"""
String that runs code: python GeneralizedLaunch.py

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
Calculated Force per motor : 1.30244e-10, Calculated Fuel Consumption per motor : 3.00233e-14
Calculated Force           : 2.17073e+06, Calculated Fuel Consumption      : 5.00388e+02

Generalized Position in solar system frame at t = 0: [x : 2.814 AU, y : 6.85e-05 AU]
Generalized Velocity in solar system frame at t = 0: [x : 2.479 AU/Y, y : 5.854 AU/Y]
Generalized Planet Position: [2.81419848 0.        ]
Generalized Planet Position: [0.       5.703728]
Generalized Position: [4.94033880e-05 1.22667618e-06]
Generalized Velocity: [2.4786584  0.15038049]

Specialized Position in solar system frame at t = 0: [x : 2.814 AU, y : 6.85e-05 AU]
Specialized Velocity in solar system frame at t = 0: [x : 2.479 AU/Y, y : 5.854 AU/Y]
Specialized Planet Position: (2.814198478717244, 0.0)
Specialized Planet Velocity: (0.0, 5.703728002035204)
Specialized Position: [4.94033880e-05 1.22667618e-06]
Specialized Velocity: [2.4786584  0.15038049]

Rocket was moved up by 1.49664e-05 m to stand on planet surface.
New launch parameters set.
Launch completed, reached escape velocity in 372.182 s.
Your spacecraft position was satisfyingly calculated. Well done!
*** Achievement unlocked: No free launch! ***

"""