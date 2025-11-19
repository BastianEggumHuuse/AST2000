# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune

# Regular imports
import numpy             as np
import matplotlib.pyplot as plt
from numba import njit
import pickle as pkl 

import Analytiskeløsninger as soloutions

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission


class LandingSimulation:

    def __init__(self,mission,t_0,r_0,v_0,dt = 10e-3):

        self.mass            = mission.lander_mass
        self.area            = mission.lander_area
        self.planet_radius   = mission.system.radii[1]*1000
        self.planet_mass     = mission.system.masses[1] * const.m_sun
        self.planet_rotation = (2*np.pi) / (mission.system.rotational_periods[1]*24*60*60)
        self.z_hat = np.array([0,0,1])

        self.initialize_density(r_0)

        self.R = np.array([r_0])
        self.V = np.array([v_0])
        self.A = np.array([self.compute_acceleration(r_0,v_0)])

        self.t  = t_0
        self.dt = dt

        self.crashed = False

    def initialize_density(self,r_0):
        r = np.linalg.norm(r_0)
        R = np.linspace(self.planet_radius,r,100000)

        T, rho, j = soloutions.TandRho(R)

        self.r_limit = R[j]
        self.T_limit = T[j]

    def drag_velocity(self,r,v):
        return v + (((r[0]**2 + r[1]**2)**0.5) * self.planet_rotation) * (np.cross(r,self.z_hat)/np.linalg.norm(r))


    def compute_acceleration(self,r,v):

        # Computing drag acceleration
        v_d = self.drag_velocity(r,v)
        v_d_hat = v_d / np.linalg.norm(v_d)
        a_d = ((0.5 * soloutions.Rho(np.linalg.norm(r),self.r_limit,self.T_limit) * self.area * np.linalg.norm(v_d)**2)/self.mass) * -v_d_hat

        # Computing gravitational acceleration
        a_g = ((const.G * self.planet_mass)/np.linalg.norm(r)**2) * (-r/np.linalg.norm(r))

        a = a_d + a_g
        return a
    
    def time_step(self,R,V,A,i):

        A[i+1] = self.compute_acceleration(R[i],V[i])
        V[i+1] = V[i] + A[i+1] * self.dt
        R[i+1] = R[i] + V[i+1] * self.dt

    def fall(self,t):

        N = int(t/self.dt)
        
        R = np.zeros((N,3))
        V = np.zeros((N,3))
        A = np.zeros((N,3))

        R[0] = self.R[-1]
        V[0] = self.V[-1]
        A[0] = self.A[-1]

        for n in range(0,N-1):
            
            # Exit clause
            if(np.linalg.norm(R[n]) < self.planet_radius and self.crashed == False):
                self.crash(f"Lander has hit the ground with velocity {np.linalg.norm(V[n])}")

            if(self.crashed):
                R[n+1] = np.array([0,0,self.planet_radius])
                V[n+1] = np.zeros(3)
                A[n+1] = np.zeros(3)

            else:
                self.time_step(R,V,A,n)

        # Updating data
        self.t += t
        self.R = np.concatenate((self.R,R[1:]))
        self.V = np.concatenate((self.V,V[1:]))
        self.A = np.concatenate((self.A,A[1:]))

    def crash(self,message):
        print(message)
        print("     _.-^^---....,,-- \n _--                  --_\n<                        >)\n|                         |\n \\._                   _./\n    ```--. . , ; .--'''\n          | |   |\n       .-=||  | |=-.\n       `-=#$%&%$#=-'\n          | ;  :|\n _____.,-#%&$@%#&#~,._____ ")
        self.crashed = True

if __name__ == "__main__":

    #Unpickling launch
    with open("Mission.pkl", 'rb') as file:
        mission = pkl.load(file)

    with open("Landing.pkl", 'rb') as file:
        landing = pkl.load(file)

    t_0,r_0,v_0 = landing.orient()
    v_0 = np.array([0,0,0])

    LandingSim = LandingSimulation(
        mission,
        t_0 = t_0,
        r_0 = r_0,
        v_0 = v_0,
        )
    
    LandingSim.fall(1000 * 4)

    r_hats = LandingSim.R / np.linalg.norm(LandingSim.R,keepdims=True,axis = 1)
    a_rs   = np.sum(LandingSim.A * r_hats,keepdims=True,axis = 1)

    ax = plt.axes()
    ax.plot(np.linspace(t_0,LandingSim.t,len(LandingSim.R)),np.linalg.norm(LandingSim.R,keepdims=True,axis = 1) - mission.system.radii[1]*1000)
    plt.show()
    ax = plt.axes()
    ax.plot(LandingSim.R[:,0],LandingSim.R[:,1])
    plt.axis("equal")
    plt.show()