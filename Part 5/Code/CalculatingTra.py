# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune

# Regular imports
import numpy             as np
import matplotlib.pyplot as plt
from numba import njit
import pickle as pkl 
import sys

from GeneralizedLaunch import NumericalOrbitFunction
# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils

from ast2000tools.space_mission import SpaceMission

FilePath = "NumericalOrbitData.npz"
PlanetPositionFunction = NumericalOrbitFunction(FilePath)

npz = np.load(FilePath)
config       = npz["config"]
TotalTime    = float(config[0])
NumSteps     = int(config[2])
OrbitTimes   = npz["OrbitTimes"]
r = npz["r"]
Info = (TotalTime, NumSteps, OrbitTimes, r)

@njit
def FindR(t,p, info):
    
    """
    Method that returns the position of a given planet along the x and y axes at a given time.

    t       : float        | the desired point in time
    p       : int          | the desired planet index

    returns : Array(float) | the position of the given planet at the given time
    """
    
    TotalTime    = info[0]
    NumSteps     = info[1]
    OrbitTimes   = info[2]
    r = info[3]
    

    # Setting total time, delta time, and number of time steps, from the read file


        # Setting r from read file
    
    # Wrapping the t-value
    # if t is less than zero, we make it wrap around to the end of the simulation
    # This stops index-issues.
    if(t < 0):
        t = OrbitTimes[p] - t
    
    # Finding the index of the given time
    # This deserves an explanation. Since the positions are stored in an array with a length of NumSteps,
    # we can't just insert t into this array to get the value (since t is a floating number)
    # t/self.TotalTime gives us the percentage of the simulation the time t is at.
    # (if t/self.TotalTime = 0.5, t is halfway through the simulation).
    # We multiply this number with the total amount of steps, to get the closes time index to our current time.
    # We then floor that index (round down) and turn it into an integer.
    Index = int(np.floor((t/TotalTime)*NumSteps))

    # Finding x and y positions at this index
    x = (r[0][p][Index])
    y = (r[1][p][Index])
    # Returning vector
    return(np.array([x,y]))

AU = const.AU
G_sol = const.G_sol

@njit
def Lerp(R_0,R_1, I):
    dR_x = R_1[0] - R_0[0]
    dR_y = R_1[1] - R_0[1]
    
    R_x = R_0[0] + dR_x * I
    R_y = R_0[1] + dR_y * I

    return R_x,R_y

@njit
def GravitationalAks(R,K,N_k,dt, M,T_0, info):
    
    if(K == 1):
        print(R[K])

    a_x = -G_sol *M[-1]*R[K][0]/ ((R[K][0])**2 + R[K][1]**2)**(3/2)
    a_y = -G_sol *M[-1]*R[K][1]/ ((R[K][0])**2 + R[K][1]**2)**(3/2)

    if(K == 1):
        print(a_x,a_y)

    for j in range(len(M)-1):
        
        R_0 = FindR(T_0 + K*dt, j, info)
                
        R_1 = FindR(T_0 + (K+N_k-1)*dt,j, info)

        R_p_x, R_p_y = Lerp(R_0, R_1, (K%N_k)/N_k) 

        r_x = R[K][0] - R_p_x
        r_y = R[K][1] - R_p_y

        gamma = -G_sol * M[j]/((r_x**2 + r_y**2)**(3/2))

        a_x += r_x * gamma
        a_y += r_y * gamma

    return a_x, a_y

@njit
def GravitationalAcceleration(R,M,dt,N_k,k,t_0,Info):
    
    # Calculating the acceleration from the sun
    r = R[k]
    r_len = (r[0]**2 + r[1]**2)**0.5
    r_hat_x = r[0]/r_len
    r_hat_y = r[1]/r_len
    
    a = -((G_sol * M[-1])/(r_len**2))
    a_x = a * r_hat_x
    a_y = a * r_hat_y

    t = t_0 + k*dt
    t_1 = t_0 + k*dt + (N_k-1) * dt
    interpolator = ((k)%N_k)/N_k

    for i in range(len(M) - 1):

        r_p_0 = FindR(t,i,Info)
        r_p_1 = FindR(t_1,i,Info)
        r_p = Lerp(r_p_0,r_p_1,interpolator)

        r_x = R[k][0] - r_p[0]
        r_y = R[k][1] - r_p[1]

        r = np.zeros(2); r[0] = r_x; r[1] = r_y
        r_len = (r[0]**2 + r[1]**2)**0.5
        r_hat_x = r[0]/r_len
        r_hat_y = r[1]/r_len

        a = -((G_sol * M[i])/(r_len**2))
        a_x += a * r_hat_x
        a_y += a * r_hat_y

    return a_x,a_y

@njit
def timestep(R,v,a,dt, N_k, K, M, T_0, info):
    R[K+1][0] = R[K][0] + v[K][0] * dt + 0.5*a[K][0]*dt**2
    R[K+1][1] = R[K][1] + v[K][1] * dt + 0.5*a[K][1]*dt**2

    #a[K+1] = GravitationalAks(R,K+1,N_k,dt, M, T_0, info)
    a[K+1] = GravitationalAcceleration(R,M,dt,N_k,K+1, T_0, info)

    if(K == 0):

        print(a[K])
        print(a[K+1])

    v[K+1][0] = v[K][0] + 0.5*(a[K][0] + a[K+1][0])*dt
    v[K+1][1] = v[K][1] + 0.5*(a[K][1] + a[K+1][1])*dt

#@njit
def Main(R_0, v_0, dt,dT, N, M,T_0, info):
    
    N_k = np.floor(dT/dt)
    R = np.zeros((N+1,2))
    R[0] = R_0

    v = np.zeros((N+1,2))
    v[0] = v_0

    a = np.zeros((N+1,2))
    #a[0] = np.array(GravitationalAks(R,0,N_k,dt, M,T_0, info))
    #a[0] = np.array(GravitationalAcceleration(R,M,dt,N_k,0, T_0, info))
    

    r = R[0] - FindR(T_0,0,info)
    r_hat = r/np.linalg.norm(r)
    a[0] = -((G_sol*M[0])/np.linalg.norm(r)) * r_hat

    r = R[0]
    r_hat = r/np.linalg.norm(r)
    a[0] += -((G_sol*M[-1])/np.linalg.norm(r)) * r_hat
    


    for K in range(N):
        timestep(R,v,a,dt,N_k, K, M,T_0, info)
    T = dt*N
    return R,v,a, T

if __name__ == '__main__':

    with open("Mission.pkl", 'rb') as file:
        mission = pkl.load(file)

    R_0 = mission._position_after_launch
    v_0 = mission._velocity_after_launch
    
    M = np.zeros(mission.system._number_of_planets + 1)
    M[:-1] = mission.system.masses
    M[-1] = mission.system.star_mass
    
    dT = PlanetPositionFunction.dt
    dt = 1/1000000
    T_0 = mission.time_after_launch# + dT

    N = 3140000
    
    R, v, a, T = Main(R_0, v_0, dt,dT, N, M ,T_0, Info)
   
    
    plt.plot(r[0][0][int(T_0/dT):int(T_0/dT)+int(N*dt/dT)],r[1][0][int(T_0/dT):int(T_0/dT)+int(N*dt/dT)] )
    R = R.T
    
    plt.plot(R[0], R[1])

    # Rocket pos
    plt.plot(R[0][:1000],R[1][:1000])

    r_p_1 = PlanetPositionFunction.range(0,10)
    plt.plot(r_p_1[0][1],r_p_1[1][1])
    plt.show()
    

