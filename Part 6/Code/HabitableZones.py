# Bruker ikke kodemal!!!
# Skrevet av Bastian Eggum Huuse og Bendik Thune

import matplotlib.pyplot as plt
import numpy as np

import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

from GeneralizedLaunch import NumericalOrbitFunction

class HabitableZones:

    def __init__(self,mission):

        """
        Class that calculates the surface temperatures for all the planets.

        mission : SpaceMission | Instance of spacemission

        returns : self
        """

        self.mission = mission
        self.system  = mission.system

        # We want to work in SI units this time around!
        self.r_star  = self.system.star_radius * 1000
        self.T_star  = self.system.star_temperature
        self.r_planets = self.system.radii * 1000

        # Getting Orbit filename and initializing position function
        self.FileName = "NumericalOrbitData.npz"
        self.R_planets = NumericalOrbitFunction(self.FileName)

        # Creating empty array which we will fill with temperatures
        self.Temps = np.zeros((self.R_planets.NumSteps,self.system.number_of_planets))

        # Getting variables from the positional function
        self.TotalTime    = self.R_planets.TotalTime
        self.dt           = self.R_planets.dt
        self.NumSteps     = self.R_planets.NumSteps

    def Loop(self):

        """
        Method that finds the temperatures for all the planets over the time self.TotalTime

        returns : Array(float) | The temperatures for each planet at every timestep
        """

        # Getting the positions of all planets for the entire duration of self.TotalTime
        # The syntax here is kind of crazy: 
        # Firstly we get the entire list of all positions. This array has dimentions (2,7,self.NumSteps), or (Axes, NumPlanets, NumSteps)
        # We transpose this array, giving us dimentions (self.NumSteps,7,2), and then take the lengths of every position.
        # This gives us dimentions (self.NumSteps,7,1), where the last dimention is the norm of the corresponding position
        # We then transpose this again, giving us our original dimentions with the first reduced to 1 : (1,7,self.NumSteps)
        # Finally we multiply by const.AU to get this in meters (the function returns in AU)
        PosRange = (np.linalg.norm((self.R_planets.range(0,self.TotalTime)).T,axis = 2, keepdims= True)).T[0] * const.AU

        # We then generate a list with dimentions (self.NumSteps,7), containing the radii of the planets repeated NumSteps times.
        # We then turn this into an array, and transpose it, giving us dimentions (7,self.NumSteps)
        # Because of this transformation, we can now do math with the above array since they have the same dimentions!
        Radii = np.array([self.r_planets] * self.NumSteps).T

        # Calculating the surface temperature for every planet, and every timestep all at once.
        # Everything that happens above is simply to let us do this step vectorized :)
        self.Temps = (self.r_star/(PosRange * 2 ))**(1/2) * self.T_star

        return(self.Temps)



    
if __name__ == "__main__":

    # Initializing AST
    Seed = utils.get_seed('bmthune')
    mission = SpaceMission(Seed)

    Zones = HabitableZones(mission)
    Temperatures = Zones.Loop()

    fig, ax = plt.subplots()

    plt.axhline(260,color = "k")
    plt.axhline(390,color = "k")
    plt.text(121,256,"260 K",color = "k")
    plt.text(121,386,"390 K",color = "k")
    plt.fill([0,120,120,0],[260,260,390,390],color = "wheat")

    i = 1
    for T in Temperatures:
        ax.plot(np.linspace(0,Zones.TotalTime,Zones.NumSteps),T,label = f"Planet {i}")
        i +=1

    plt.legend(loc = 'upper right')
    plt.xlabel("Tid [Y]")
    plt.ylabel("Temperatur [K]")
    plt.title("Temperatur for planetene i solsystemet over tid")

    plt.xlim(0,120)

    plt.show()

"""
String that runs code: python HabitableZones.py

This code has no text output, but does output an image. This image has been turned in with the name "HabitableZonesPlot.png".
"""
