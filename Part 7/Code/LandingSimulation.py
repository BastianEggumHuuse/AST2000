# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune

# Regular imports
import numpy             as np
import matplotlib.pyplot as plt
from numba import njit
import pickle as pkl 

import Analytiskeløsninger as soloutions
from LandingPreparation import CoordinateAtTime,ComputeAngle

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission


class LandingSimulation:

    def __init__(self,mission,t_0,r_0,v_0,dt = 10e-3):

        """

        LandingSimulation is a class written to imitate the way the ast2000Tools has
        written it's classes. We initialize an instance, and then have multiple choises on 
        how the instance will operate, for example falling for t seconds, or opening the parachute
        
        parameters:
        mission : SpaceMission | Instance of the SpaceMission class
        t_0     : float        | Initial time in seconds
        r_0     : Array(float) | Initial positional vector
        v_0     : Array(float) | Initial velocity vector
        dt      : float        | size of time steps

        returns:
        self
        """

        print(f"\nInitializing landing simulation at t = {t_0}\n")

        # Collecting information about the system
        self.mass            = mission.lander_mass
        self.lander_area     = mission.lander_area
        self.lander_limit    = 1e7 #Pa
        self.parachute_area  = 13 # m^2
        self.parachute_limit = 250000 # N
        self.planet_radius   = mission.system.radii[1]*1000
        self.planet_mass     = mission.system.masses[1] * const.m_sun
        self.planet_rotation = (2*np.pi) / (mission.system.rotational_periods[1]*24*60*60)
        self.z_hat           = np.array([0,0,1])

        # Storing lander states
        self.parachute_open   = False
        self.parachute_broken = False
        self.landed           = False

        # Initializing density
        self._initialize_density(r_0)

        # Initializing arrays
        self.R = np.array([r_0])
        self.V = np.array([v_0])
        self.A = np.array([self._compute_acceleration(r_0,v_0, 0)])

        # Initializing time
        # A lot of different t values are stored here that aren't mentioned in the flowchart.
        # We realized that we wanted to also keep track of how long the simulation has been running,
        # when we crash/land, etc.
        self.t_0     = t_0
        self.t       = self.t_0
        self.sim_t   = 0
        self.final_t = t_0
        self.dt      = dt

    def _initialize_density(self,r_0):

        """
        This method was not detailed in the flowchart, as we realized we needed it after writing said flowchart.
        We have to find the value r_limit where the atmosphere becomes adiabatic

        parameters:
        r_0 : Array(float) | the initial positional vector
        """

        r = np.linalg.norm(r_0)
        R = np.linspace(self.planet_radius,r,100000)

        # Finding rho (from part 6)
        T, rho, j = soloutions.TandRho(R)

        self.r_limit = R[j]
        self.T_limit = T[j]
        print("Density Initialized.\n")

    def _drag_velocity(self,r,v):
        """
        Method that computes the v_drag of the lander

        parameters:
        r : Array(float) | current position
        v : Array(float) | current velocity

        returns:
        Array(float) | v_drag
        """

        return v - (((r[0]**2 + r[1]**2)**0.5) * self.planet_rotation) * (np.cross(self.z_hat,r)/np.linalg.norm(r))


    def _compute_acceleration(self,r,v,n):

        """
        Method that computes the total acceleration of the lander at a point in time

        parameters:
        r : Array(float) | current position
        v : Array(float) | current velocity

        returns:
        Array(float) | total acceleration of the lander
        """

        # Computing drag velocity
        v_d = self._drag_velocity(r,v)
        v_d_hat = v_d / np.linalg.norm(v_d)

        # Computing drag acceleration
        if(self.parachute_open):

            # If the parachute is open, we use a larger area
            a_d = ((0.5 * soloutions.Rho(np.linalg.norm(r),self.r_limit,self.T_limit) * self.parachute_area * np.linalg.norm(v_d)**2)/self.mass) * -v_d_hat
            
            # Checking if the parachute breaks
            if(np.linalg.norm(a_d) * self.mass > self.parachute_limit):
                print(f"Drag force is {np.linalg.norm(a_d) * self.mass} which exceeds {self.parachute_limit}!")
                print("The parachute is now broken.")
                self.parachute_broken = True
                self.parachute_open   = False
            
        else:        
            # If the parachute isn't open, we use the normal, smaller area   
            a_d = ((0.5 * soloutions.Rho(np.linalg.norm(r),self.r_limit,self.T_limit) * self.lander_area * np.linalg.norm(v_d)**2)/self.mass) * -v_d_hat
        
        # Calculating the drag pressure the lander experiences
        P_d = ((0.5 * soloutions.Rho(np.linalg.norm(r),self.r_limit,self.T_limit) * np.linalg.norm(v_d)**2))
        
        # Checking if the lander burns
        if P_d > self.lander_limit:
            self.final_t = self.t + n*self.dt
            self.burn()

        # Computing gravitational acceleration
        a_g = ((const.G * self.planet_mass)/np.linalg.norm(r)**2) * (-r/np.linalg.norm(r))

        # Computing total acceleration
        a = a_d + a_g
        return a
    
    def _time_step(self,R,V,A,i):

        """
        Updating arrays for a given time step i

        parameters:
        R : Array(float) | Array of positions
        V : Array(float) | Array of velocities
        A : Array(float) | Array of accelerations
        i : int          | index of current time step
        """

        A[i+1] = self._compute_acceleration(R[i],V[i],i)
        V[i+1] = V[i] + A[i+1] * self.dt
        R[i+1] = R[i] + V[i+1] * self.dt

    def fall(self,t):

        """
        The main method we use to interact with the class.
        This method simulates a fall through the atmosphere for t seconds.

        parameters:
        t : float | how long to fall for
        """

        print(f"Falling for {t} seconds.")

        # Computing num time steps
        N = int(t/self.dt)
        
        # Initializing temporary arrays
        R = np.zeros((N,3))
        V = np.zeros((N,3))
        A = np.zeros((N,3))

        # Initializing array values
        R[0] = self.R[-1]
        V[0] = self.V[-1]
        A[0] = self.A[-1]

        # Looping over N time steps
        for n in range(0,N-1):
            
            # Exit clause
            # If this test succedes, we've either landed or crashed
            if(np.linalg.norm(R[n]) <= self.planet_radius and self.landed == False):
                self.final_t = self.t + n*self.dt
                
                # Computing our velocity in relation to the ground
                ground_velocity = V[n] - (((R[n][0]**2 + R[n][1]**2)**0.5) * self.planet_rotation) * (np.cross(self.z_hat,R[n])/np.linalg.norm(R[n]))

                # Checking if we've landed or crashed
                if(np.linalg.norm(ground_velocity) < 3):
                    self.land(f"Lander has hit the ground with velocity {np.linalg.norm(ground_velocity)}.")
                else:
                    self.crash(f"Lander has hit the ground with velocity {np.linalg.norm(ground_velocity)}.")

            # If we've landed (or crashed) the arrays stay constant
            # this is for plotting reasons
            if(self.landed):
                R[n+1] = R[n]
                V[n+1] = V[n]
                A[n+1] = A[n]
            # If we're still falling, compute stuff :)
            else:
                self._time_step(R,V,A,n)

        # Updating data
        self.t += t
        self.sim_t += t
        # Adding the entirity of the temp arrays to our main arrays
        self.R = np.concatenate((self.R,R[1:]))
        self.V = np.concatenate((self.V,V[1:]))
        self.A = np.concatenate((self.A,A[1:]))

    def open_parachute(self):

        """
        This method simply sets parachute_open to true
        """

        if(self.parachute_broken == False):
            print("Opening parachute.")
            self.parachute_open = True

    def burn(self):

        """
        Method to signify the lander burning
        """

        print(f"You burned to a crisp t = {self.final_t}, sim_time {self.final_t - self.t_0}")
        self.landed = True

    def land(self,message):
        """
        Method to signify the lander landing.

        parameters:
        message : string | The message to print when the lander lands
        """

        print(f"A landing has occured at t = {self.final_t}, sim_time {self.final_t - self.t_0}")
        print(message)
        print("_--^*# Landed Succesfully #*^--_")
        print("                     `. ___                                  \n                    __,' __`.                _..----....____ \n        __...--'``;.   ,.   ;``--..__     .'    ,-._    _.-'\n  _..-''-------'   `'   `'   `'     O ``-''._   (,;') _,'    \n,'________________                          \\`-._`-','       \n `._              ```````````------...___   '-.._'-:         \n    ```--.._      ,.                     ````--...__\\-.      \n            `.--. `-`                       ____    |  |`    \n              `. `.                       ,'`````.  ;  ;`    \n                `._`.        __________   `.      \\'__/`     \n                   `-:._____/______/___/____`.     \\  `      \n                               |       `._    `.    \\        \n                               `._________`-.   `.   `.___   \n                                             SSt  `------'`' ")
        self.landed = True

    def crash(self,message):
        """
        Method to signify the lander crashing.

        parameters:
        message : string | The message to print when the lander crashes
        """

        print(f"A crash has occured at t = {self.final_t}, sim_time {self.final_t - self.t_0}")
        print(message)
        print("     _.-^^---....,,-- \n _--                  --_\n<                        >)\n|                         |\n \\._                   _./\n    ```--. . , ; .--'''\n          | |   |\n       .-=||  | |=-.\n       `-=#$%&%$#=-'\n          | ;  :|\n _____.,-#%&$@%#&#~,._____ ")
        self.landed = True

if __name__ == "__main__":

    #Unpickling launch
    with open("Mission.pkl", 'rb') as file:
        mission = pkl.load(file)

    #Unpickling landing
    with open("Landing.pkl", 'rb') as file:
        landing = pkl.load(file)

    t_1,r_1,v_1 = landing.orient()
    landing.fall(1200)
    t_0,r_0,v_0 = landing.orient()
    v_0 = np.array([0,0,-180]) -861*(np.cross(r_0,np.array([0,0,1])))/np.linalg.norm(r_0)


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
    LandingSim.fall(4000)

    # Info time :)

    # Turning start position into spherical coordinates
    start_position_r = np.linalg.norm(r_0)
    start_position_phi = ComputeAngle(r_0)
    start_position_theta = np.arccos(r_0[2]/np.linalg.norm(r_0))
    position_t_0 = np.array([start_position_r,start_position_phi,start_position_theta])
    # Turning end position into spherical coordinates
    lander_position_r     = np.linalg.norm(LandingSim.R[-1])
    lander_position_phi   = ComputeAngle(LandingSim.R[-1])
    lander_position_theta = np.arccos(LandingSim.R[-1][2]/np.linalg.norm(LandingSim.R[-1]))
    lander_position_t     = np.array([lander_position_r,lander_position_phi,lander_position_theta])

    destination_position_0 = np.array([2304594.3015970597,4.8107257594915245,1.6083942634485817])#[2304594.3015970597,3.5822217279404853,1.608027384048026])
    destination_position_t = CoordinateAtTime(destination_position_0,LandingSim.final_t - LandingSim.t_0,LandingSim.planet_rotation)
    
    print(np.rad2deg(ComputeAngle(r_1)))
    print(np.rad2deg(start_position_phi))
    print('Ship(t = 0) : ', position_t_0)
    print('Ship(t = 0) : ', destination_position_0)
    print("Lander      : ", lander_position_t)
    print("Destination : ", destination_position_t)

    # Plotting
    ax = plt.axes()
    ax.plot(np.linspace(0, LandingSim.sim_t,len(LandingSim.R)),np.linalg.norm(LandingSim.V,keepdims=True,axis = 1))
    plt.axvline(LandingSim.final_t-t_0)
    plt.show()

    ax = plt.axes()
    ax.plot(np.linspace(0,LandingSim.sim_t,len(LandingSim.R)),np.linalg.norm(LandingSim.R,keepdims=True,axis = 1) - mission.system.radii[1]*1000)
    plt.xlabel("Tid [s]")
    plt.ylabel("Distanse fra overflaten [m]")
    plt.axvline(LandingSim.final_t-t_0,color = "green")
    plt.axvline(800,color = "red")
    plt.show()

    ax = plt.axes()
    ax.plot(LandingSim.R[:,0],LandingSim.R[:,1])
    Planet = plt.Circle((0,0), LandingSim.planet_radius, color = 'blue', fill = False, ls = '-')
    ax.add_patch(Planet)
    Aro = OrbitRadi = plt.Circle((0,0), LandingSim.r_limit, color = 'red', fill = False, ls = '-')
    ax.add_patch(Aro)
    plt.xlabel("Posisjon langs x-aksen [m]")
    plt.ylabel("Posisjon langs y-aksen [m]")

    plt.axis("equal")
    plt.show()

r"""

Initializing landing simulation at t = 82980.0

c:\Users\tmthu\OneDrive\Documents\Bendik\UiO\AST2200\Prosjekt\Git3\AST2000\Part 7\Code\Analytiskeløsninger.py:51: RuntimeWarning: invalid value encountered in power
  rho = (2/5*(-a*g/gamma * (mu*const.m_p/const.k_B)**gamma * r + C))**(5/2)
Density Initialized.

Falling for 800 seconds.
Opening parachute.
Falling for 4000 seconds.
A landing has occured at t = 86970.04, sim_time 3990.0399999999936
Lander has hit the ground with velocity 2.9784796244091205.
_--^*# Landed Succesfully #*^--_
                     `. ___
                    __,' __`.                _..----....____
        __...--'``;.   ,.   ;``--..__     .'    ,-._    _.-'
  _..-''-------'   `'   `'   `'     O ``-''._   (,;') _,'
,'________________                          \`-._`-','
 `._              ```````````------...___   '-.._'-:
    ```--.._      ,.                     ````--...__\-.
            `.--. `-`                       ____    |  |`
              `. `.                       ,'`````.  ;  ;`
                `._`.        __________   `.      \'__/`
                   `-:._____/______/___/____`.     \  `
                               |       `._    `.    \
                               `._________`-.   `.   `.___
                                             SSt  `------'`'
Ship(t = 0) :  [2.68468998e+06 4.65440989e+00 1.57079633e+00]
Ship(t = 0) :  [2.30459430e+06 4.81072576e+00 1.60839426e+00]
Lander      :  [2.30459427e+06 4.99915807e+00 1.60835421e+00]
Destination :  [2.30459430e+06 4.99924734e+00 1.60839426e+00]
"""
