# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune

# Regular imports
import numpy             as np
import matplotlib.pyplot as plt

# Orbit imports
from OrbitPlotNumerical import NumericalOrbitFunction

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

class KeplerArea:

    def __init__(self,Filename,PlanetIndex):

        self.Filename = Filename
        self.PlanetIndex = PlanetIndex

        self.r = NumericalOrbitFunction(self.Filename)
        self.dt = self.r.dt

    def __call__(self, t0 = 0, t1 = 1):
        
        A = 0

        t = t0
        while t < t1:
            r0 = self.r(t,self.PlanetIndex)
            t += self.dt
            r1 = self.r(t,self.PlanetIndex)

            # COMMENT CROSS PRODUCT
            dA = 0.5 * abs(np.cross(r0,r1))
            A += dA

        return A
    
    
if __name__ == "__main__":

    # Initializing AST
    Seed = utils.get_seed('bmthune')
    mission = SpaceMission(Seed)
    system = mission.system

    R0 = system.initial_positions
    V0 = system.initial_velocities

    # Initializing Area Function
    AreaFunction = KeplerArea("NumericalOrbitData.npz",0)
    TotalTime = AreaFunction.r.OrbitTimes[0]

    # Finding radius at aphelion and perihelion COMMENT THIS!!!!!!!!!
    t_aph = system.aphelion_angles[0]/(2*np.pi) * AreaFunction.r.OrbitTimes[0]
    t_per = (system.aphelion_angles[0] + np.pi)/(2*np.pi) * AreaFunction.r.OrbitTimes[0]

    t_interval = 0.05
    
    # Finding the positions at 
    r = AreaFunction.r.range(AreaFunction.dt,TotalTime)   
    r_aph = AreaFunction.r.range(t_aph - t_interval,t_aph + t_interval)
    r_per = AreaFunction.r.range(t_per - t_interval,t_per + t_interval)

    r_aph_polygon = np.concatenate((np.zeros((2,7,1)),r_aph,np.zeros((2,7,1))),2)
    r_per_polygon = np.concatenate((np.zeros((2,7,1)),r_per,np.zeros((2,7,1))),2)

    # Initializing plotting
    fig, ax = plt.subplots()

    ax.plot(r[0][0],r[1][0],color = "gold")
    ax.fill(r_aph_polygon[0][0],r_aph_polygon[1][0],color = "royalblue")
    ax.fill(r_per_polygon[0][0],r_per_polygon[1][0],color = "red")

    A_aph = AreaFunction(t_aph - t_interval,t_aph + t_interval)
    A_per = AreaFunction(t_per - t_interval,t_per + t_interval)

    print(f"Aphelion area : {A_aph:.5f} AU^2 | Perihelion area : {A_per:.5f} AU^2 | Ratio = {A_aph/A_per:.5f}")
    

    # Making axes equal, and showing plot
    plt.axis('equal')
    plt.show()
    