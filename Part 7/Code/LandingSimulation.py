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

        print(f"\nInitializing landing simulation at t = {t_0}\n")

        self.mass            = mission.lander_mass
        self.lander_area     = mission.lander_area
        self.lander_limit    = 1e7 #Pa
        self.parachute_area  = 13 # m^2
        self.parachute_limit = 250000 # N
        self.planet_radius   = mission.system.radii[1]*1000
        self.planet_mass     = mission.system.masses[1] * const.m_sun
        self.planet_rotation = (2*np.pi) / (mission.system.rotational_periods[1]*24*60*60)
        self.z_hat           = np.array([0,0,1])

        self.parachute_open   = False
        self.parachute_broken = False
        self.landed           = False

        self._initialize_density(r_0)

        self.R = np.array([r_0])
        self.V = np.array([v_0])
        self.A = np.array([self._compute_acceleration(r_0,v_0)])

        self.t_0     = t_0
        self.t       = self.t_0
        self.sim_t   = 0
        self.final_t = 0
        self.dt      = dt

    def _initialize_density(self,r_0):
        r = np.linalg.norm(r_0)
        R = np.linspace(self.planet_radius,r,100000)

        T, rho, j = soloutions.TandRho(R)

        self.r_limit = R[j]
        self.T_limit = T[j]
        print("Density Initialized.\n")

    def _drag_velocity(self,r,v):
        return v - (((r[0]**2 + r[1]**2)**0.5) * self.planet_rotation) * (np.cross(r,self.z_hat)/np.linalg.norm(r))


    def _compute_acceleration(self,r,v):

        # Computing drag acceleration
        v_d = self._drag_velocity(r,v)
        v_d_hat = v_d / np.linalg.norm(v_d)

        # Computing drag acceleration
        if(self.parachute_open):
            a_d = ((0.5 * soloutions.Rho(np.linalg.norm(r),self.r_limit,self.T_limit) * self.parachute_area * np.linalg.norm(v_d)**2)/self.mass) * -v_d_hat
            
            if(np.linalg.norm(a_d) * self.mass > self.parachute_limit):
                print(f"Drag force is {np.linalg.norm(a_d) * self.mass} which exceeds {self.parachute_limit}!")
                print("The parachute is now broken.")
                self.parachute_broken = True
                self.parachute_open   = False
            
        else:           
            a_d = ((0.5 * soloutions.Rho(np.linalg.norm(r),self.r_limit,self.T_limit) * self.lander_area * np.linalg.norm(v_d)**2)/self.mass) * -v_d_hat
        
        P_d = ((0.5 * soloutions.Rho(np.linalg.norm(r),self.r_limit,self.T_limit) * np.linalg.norm(v_d)**2)/self.mass)
        
        if P_d > self.lander_limit:
            self.final_t = self.t + n*self.dt
            self.burn()
        # Computing gravitational acceleration
        a_g = ((const.G * self.planet_mass)/np.linalg.norm(r)**2) * (-r/np.linalg.norm(r))

        a = a_d + a_g
        return a
    
    def _time_step(self,R,V,A,i):

        A[i+1] = self._compute_acceleration(R[i],V[i])
        V[i+1] = V[i] + A[i+1] * self.dt
        R[i+1] = R[i] + V[i+1] * self.dt

    def fall(self,t):

        print(f"Falling for {t} seconds.")

        N = int(t/self.dt)
        
        R = np.zeros((N,3))
        V = np.zeros((N,3))
        A = np.zeros((N,3))

        R[0] = self.R[-1]
        V[0] = self.V[-1]
        A[0] = self.A[-1]

        for n in range(0,N-1):
            
            # Exit clause
            if(np.linalg.norm(R[n]) <= self.planet_radius and self.landed == False):
                self.final_t = self.t + n*self.dt
                Ground_velocity = V[n] - (((R[n][0]**2 + R[n][1]**2)**0.5) *self.planet_rotation) * (np.cross(R[n],self.z_hat)/np.linalg.norm(R[n]))

                if(np.linalg.norm(Ground_velocity) < 3):
                    self.land(f"Lander has hit the ground with velocity {np.linalg.norm(Ground_velocity)}.")
                else:
                    self.crash(f"Lander has hit the ground with velocity {np.linalg.norm(Ground_velocity)}.")

            
            if(self.landed):
                R[n+1] = np.array([0,0,self.planet_radius])
                V[n+1] = np.zeros(3)
                A[n+1] = np.zeros(3)
            
            else:
                self._time_step(R,V,A,n)

        # Updating data
        self.t += t
        self.sim_t += t
        self.R = np.concatenate((self.R,R[1:]))
        self.V = np.concatenate((self.V,V[1:]))
        self.A = np.concatenate((self.A,A[1:]))

    def open_parachute(self):
        if(self.parachute_broken == False):
            print("Opening parachute.")
            self.parachute_open = True

    def boost(self,dv):
        self.V[-1] += dv
    def burn(self):
        print(f"You burned to a crisp t = {self.final_t}, sim_time {self.final_t - self.t_0}")
        self.landed = True

    def land(self,message):
        print(f"A landing has occured at t = {self.final_t}, sim_time {self.final_t - self.t_0}")
        print(message)
        print("_--^* Landed Succesfully *^--_")
        self.landed = True

    def crash(self,message):
        print(f"A crash has occured at t = {self.final_t}, sim_time {self.final_t - self.t_0}")
        print(message)
        print("     _.-^^---....,,-- \n _--                  --_\n<                        >)\n|                         |\n \\._                   _./\n    ```--. . , ; .--'''\n          | |   |\n       .-=||  | |=-.\n       `-=#$%&%$#=-'\n          | ;  :|\n _____.,-#%&$@%#&#~,._____ ")
        self.landed = True

if __name__ == "__main__":

    #Unpickling launch
    with open("Mission.pkl", 'rb') as file:
        mission = pkl.load(file)

    with open("Landing.pkl", 'rb') as file:
        landing = pkl.load(file)

    t_0,r_0,v_0 = landing.orient()
    v_0 = np.array([1,0,0])#(2*np.pi) / (mission.system.rotational_periods[1]*24*60*60)*(np.cross(r_0,np.array([0,0,1])))#np.array([0,0,0])
    print(np.linalg.norm(v_0))
    # Initializing simulation
    LandingSim = LandingSimulation(
        mission,
        t_0 = t_0,
        r_0 = r_0,
        v_0 = v_0,
        )
    
    # Running simulation
    LandingSim.fall(800)
    LandingSim.open_parachute()
    LandingSim.fall(3000)

    # Plotting
    ax = plt.axes()
    ax.plot(np.linspace(0,LandingSim.sim_t,len(LandingSim.R)),np.linalg.norm(LandingSim.V,keepdims=True,axis = 1))
    plt.show()
    ax = plt.axes()
    ax.plot(LandingSim.R[:,0],LandingSim.R[:,1])
    plt.axis("equal")
    plt.show()