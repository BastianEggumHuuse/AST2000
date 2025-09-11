import numpy as np
import matplotlib.pyplot as plt
# AST imports
import ast2000tools.constants as const
import ast2000tools.utils as utils
from ast2000tools.space_mission import SpaceMission


class Orbitsim:
    def __init__ (self,mission, const, TotalTime, InitialPos, InitialVel):
        #Storning AST imports
        self.mission = mission
        self.system = mission.system 
        self.const = const
        self.G = const.G_sol
        self.SM = const.m_sun

        #Time
        self.T = TotalTime
        self.dt = 1/10000

        #Initial values
        self.r0 = InitialPos.T
        self.v0 = InitialVel.T

        #Arrays
        self.N_timesteps = int(TotalTime/self.dt)
        self.r = np.zeros((self.N_timesteps, len(self.r0), 2))
        self.v = np.copy(self.r)
        self.a = np.copy(self.r)

        self.r[0] = self.r0
        self.v[0] = self.v0
        self.a[0] = ((self.G*self.SM)/self.norm(self.r[0])**2)*self.hat(self.r[0])
        self.t = np.linspace(0,self.T,self.N_timesteps)    

    def norm(self, v):
        return np.linalg.norm(v, axis=1, keepdims=True)
    def hat(self, v):
        return v/self.norm(v)

    def timestep(self,i):
        
        self.r[i] = self.r[i-1] + self.v[i-1]*self.dt + 0.5 * self.a[i-1]*self.dt**2 
        self.a[i] = ((self.G*self.SM)/self.norm(self.r[i])**2)*self.hat(self.r[i])
        self.v[i] = self.v[i-1] + 0.5*(self.a[i-1] + self.a[i]) *self.dt
    
    def loop(self):
        for i in range(1,len(self.t)):
            self.timestep(i)

        return(self.r,self.v,self.a,self.t)
    

if __name__ == "__main__":
    Seed = utils.get_seed('bmthune')

    mission = SpaceMission(Seed)
    system = mission.system

    R0 = system.initial_positions
    V0 = system.initial_velocities
    TotalTime = 20*system.rotational_periods[0]*365.25

    OrbitPlot = Orbitsim(mission,const, TotalTime, R0, V0)

    r,v,a,t = OrbitPlot.loop()
    r = r.T

    fig, ax = plt.subplots()
    # 1 = 1, 2= 2, 3= 5 , 4 = 4, 5= ?, 6= ? 7=3, 
    colours = ["red", 'darkorange', 'mediumblue', 'limegreen', 'purple', 'darkviolet', 'gold']
    #plot
    for p in range(len(OrbitPlot.r0)):
        print(p)
        ax.plot(r[0][p],r[0][p], color = colours[p])

    sol = plt.Circle((0, 0), 1, color = 'gold')
    plt.axis('equal')
    ax.add_patch(sol)
    plt.show()