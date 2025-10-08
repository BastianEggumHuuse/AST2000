import sys
sys.path.insert(0, "../../Part 2/Code")

import matplotlib.pyplot as plt
import numpy as np

import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

from OrbitPlotNumerical import NumericalOrbitFunction

class HabitableZones:

    def __init__(self,mission):

        self.mission = mission
        self.system  = mission.system

        # We want to work in SI units this time around!
        self.r_star  = self.system.star_radius * 1000
        self.T_star  = self.system.star_temperature
        self.r_planets = self.system.radii * 1000

        self.FileName = "NumericalOrbitData.npz"
        self.R_planets = NumericalOrbitFunction(self.FileName)

        self.Temps = np.zeros((self.R_planets.NumSteps,self.system.number_of_planets))

        self.TotalTime    = self.R_planets.TotalTime
        self.dt           = self.R_planets.dt
        self.NumSteps     = self.R_planets.NumSteps

    def Loop(self):

        PosRange = (np.linalg.norm((self.R_planets.range(0,self.TotalTime)).T,axis = 2, keepdims= True)).T[0] * const.AU
        Radii = np.array([self.r_planets] * self.NumSteps).T

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
