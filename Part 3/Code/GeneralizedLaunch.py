import sys
sys.path.insert(0, "../..")

import matplotlib.pyplot as plt
import numpy as np

import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

sys.path.insert(0, "../../Part 1/Code")

from SimulationRocket import SimulationRocket

sys.path.insert(0, "../../Part 2/Code")

from OrbitPlotNumerical import NumericalOrbitFunction

class GeneralizedRocket(SimulationRocket):

    def __init__(self,mission,FuelMass,SpeedBoost,NumMotors,NumParticles = 10**5,dt = 10**(-3)):
        super().__init__(FuelMass,SpeedBoost,NumMotors,NumParticles,dt)

        # Tracking some stats about the planet
        self.mission = mission
        self.system = self.mission.system

        self.FileName = "NumericalOrbitData.npz"
        self.r_planets = NumericalOrbitFunction(self.FileName)

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

        # We have a position on the y axis, radially outward from the planet, and a position on the x-axis, normally along this axis.
        # Both of these are in meters, so we must change to AU:
        self.Position *= (1/const.AU)
        self.Velocity *= ((60*60*24*365)/const.AU)

        # Now we can calculate the velocity.
        # First we get the current velocity of the planet:
        v_p = self.r_planets.GetVelocity(t_0,0)
        # Then, we swap the velocity axes around
        self.Velocity[0], self.Velocity[1] = self.Velocity[1], self.Velocity[0]
        # And finally add the velocities together:
        SolarSystemVel = v_p + np.array([self.Velocity[0] * np.cos(theta) - self.Velocity[1] * np.sin(theta), self.Velocity[0] * np.sin(theta) + self.Velocity[1] * np.cos(theta)])

        # Now, we want to get the current planet position.
        r_p = self.r_planets(t_0,0)
        # Note that we swap the x and y axes here (just like the original program)
        self.Position[0], self.Position[1] = self.Position[1], self.Position[0]
        # We rotate the local positional vector by the angle theta (Notice that this is just multiplying with the rotation matrix), and add the solar system position of the planet.
        SolarSystemPos = r_p + np.array([self.Position[0] * np.cos(theta) - self.Position[1] * np.sin(theta), self.Position[0] * np.sin(theta) + self.Position[1] * np.cos(theta)])
        # Adding the contribution of the planet's orbit speed
        SolarSystemPos += v_p * (self.t / (60 * 60 * 24 * 365))

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

    # Creating a generalized rocket, and a simulation rocket.
    GenRocket = GeneralizedRocket(mission=mission,FuelMass=Fuel,SpeedBoost=EscapeVelocity,NumMotors=NumMotors,NumParticles=Particles)
    SimRocket = SimulationRocket(FuelMass=Fuel,SpeedBoost=EscapeVelocity,NumMotors=NumMotors,NumParticles=Particles)

    # Looping both rockets
    GenRocket.TimeLoop()
    SimRocket.TimeLoop()

    Gen_Pos, Gen_Vel = GenRocket.SolarSystemPosition(0,0)
    Sim_Pos, Sim_Vel = SimRocket.StarPosition(SimRocket.Position,SimRocket.Velocity)

    print(f"Generalized Position in solar system frame at t = 0: [x : {Gen_Pos[0]:.3f} AU, y : {Gen_Pos[1]:.2e} AU]")
    print(f"Generalized Velocity in solar system frame at t = 0: [x : {Gen_Vel[0]:.3f} AU/Y, y : {Gen_Vel[1]:.3f} AU/Y]")

    print(f"Specialized Position in solar system frame at t = 0: [x : {Sim_Pos[0]:.3f} AU, y : {Sim_Pos[1]:.2e} AU]")
    print(f"Specialized Velocity in solar system frame at t = 0: [x : {Sim_Vel[0]:.3f} AU/Y, y : {Sim_Vel[1]:.3f} AU/Y]")

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
