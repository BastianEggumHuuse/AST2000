# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune

# Regular imports
import numpy             as np
import matplotlib.pyplot as plt
from numba import njit
import pickle as pkl 

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
    dR_x = -R_0[0] + R_1[0]
    dR_y = -R_0[1] + R_1[1]
    
    return R_0[0] + I* dR_x, R_0[1] + I*dR_y

@njit
def GravitationalAks(R,K,N_k,dt, M,T_0, info):
    
    a_x = -G_sol *M[-1]*R[K][0]/ (((R[K][0])**2 + R[K][1]**2)**(3/2))
    a_y = -G_sol *M[-1]*R[K][1]/ (((R[K][0])**2 + R[K][1]**2)**(3/2))

    
    for j in range(len(M)-1):
        
        R_0 = FindR(T_0 + K*dt, j, info)
                
        R_1 = FindR(T_0 + (K+N_k-1)*dt,j, info)

        R_p_x, R_p_y = Lerp(R_0, R_1, (K%N_k)/N_k) 

        r_x = -R[K][0] + R_p_x
        r_y = -R[K][1] + R_p_y

        gamma = G_sol * M[j]/((r_x**2 + r_y**2)**(3/2))

            
        a_x += r_x * gamma
        a_y += r_y * gamma

    return a_x, a_y

def GravitationAccelerations(R,K,N_k,dt, M,T_0, info):

    A = []
    Names = []

    a_x = -G_sol *M[-1]*R[K][0]/ (((R[K][0])**2 + R[K][1]**2)**(3/2))
    a_y = -G_sol *M[-1]*R[K][1]/ (((R[K][0])**2 + R[K][1]**2)**(3/2))
    a_sol = np.array([a_x,a_y])
    Names.append("Sol")
    A.append(np.linalg.norm(a_sol))

    for j in range(len(M)-1):
        
        R_0 = FindR(T_0 + K*dt, j, info)
                
        R_1 = FindR(T_0 + (K+N_k-1)*dt,j, info)

        R_p_x, R_p_y = Lerp(R_0, R_1, (K%N_k)/N_k) 

        r_x = -R[K][0] + R_p_x
        r_y = -R[K][1] + R_p_y

        gamma = G_sol * M[j]/((r_x**2 + r_y**2)**(3/2))

            
        a_x += r_x * gamma
        a_y += r_y * gamma
        a = np.array([r_x * gamma,r_y * gamma])
        A.append(np.linalg.norm(a))
        Names.append(f"Planet {j + 1}")

    return A,Names


@njit
def timestep(R,v,a,dt, N_k, K, M, T_0, info):
    R[K+1][0] = R[K][0] + v[K][0] * dt + 1/2*a[K][0]*dt**2
    R[K+1][1] = R[K][1] + v[K][1] * dt + 1/2*a[K][1]*dt**2

    a[K+1] = GravitationalAks(R,K+1,N_k,dt, M, T_0, info)
    
    v[K+1][0] = v[K][0] + 0.5*(a[K][0] + a[K+1][0])*dt
    v[K+1][1] = v[K][1] + 0.5*(a[K][1] + a[K+1][1])*dt
 
    

@njit
def Main(R_0, v_0, dt,dT, n, M,T_0, info):
    N = n+1

    
    N_k = np.floor(dT/dt)
    R = np.zeros((N,2))
    R[0] = R_0
    

    v = np.zeros((N,2))
    v[0] = v_0

    a = np.zeros((N,2))
    a[0] = np.array(GravitationalAks(R,0,N_k,dt, M,T_0, info))

    for K in range(N-1):
        timestep(R,v,a,dt,N_k, K, M,T_0, info)
        
    T = N*dt
    return R,v,a, T

if __name__ == '__main__':
    with open("Mission.pkl", 'rb') as file:
        mission = pkl.load(file)
    R_0 = mission._position_after_launch
    
    # Interpolation Tests
    print("---  Interpolation  Tests  ---")

    Interpolations = [
        [np.array([0,0]),np.array([1,0]),0.4,np.array([0.4,0])],
        [np.array([2,0]),np.array([3,0]),0.4,np.array([2.4,0])],
        [np.array([2,5]),np.array([3,4]),0.7,np.array([2.7,4.3])],
        [np.array([2313.023,5283.382]),np.array([1000.2,1392.2321]),0.5,np.array([1656.6115,3337.80705])],
        ]

    for i in Interpolations:
        print(f"Start : {i[0]}, End : {i[1]}, Interpolator : {i[2]}, Expected Value : {i[3]}, Output : {Lerp(i[0],i[1],i[2])}")

    print("")

    for i in Interpolations:
        Lerped = np.array(Lerp(i[0],i[1],i[2])) - i[0]
        Lerped = np.linalg.norm(Lerped) / np.linalg.norm(i[1] - i[0])
        Lerped = np.linalg.norm(Lerped)
        print(f"Start : {i[0]}, End : {i[1]}, Interpolator : {i[2]}, Percent of diff : {Lerped:.5f}")

    print("--- Finished Interpolation --- \n\n")

    v_0 =mission._velocity_after_launch
    M = np.zeros(mission.system._number_of_planets + 1)
    M[:-1] = mission.system.masses
    M[-1] = mission.system.star_mass
    T_0 = mission.time_after_launch
    dT = PlanetPositionFunction.dt
    dt = 1/1000000
    T_0 = mission.time_after_launch# + dT

    t = 3.14
    N = int(t / dt)#3140000 * 2
    
    R, v, a, T = Main(R_0, v_0, dt,dT, N, M ,T_0, Info)
    R_T = R.T
    ranged = PlanetPositionFunction.range(0,6)

    print("---    Gravity    Tests    ---")

    A_0,Names_0 = GravitationAccelerations(R,0,np.floor(dT/dt),dt, M, T_0, Info)
    A_1,Names_1 = GravitationAccelerations(R,int(0.5 / dt),np.floor(dT/dt),dt, M, T_0, Info)

    # Getting the first
    A_0_first = max(A_0)
    A_0_index = A_0.index(A_0_first)
    Name_0_first = Names_0[A_0_index]
    A_0.remove(A_0_first)
    Names_0.remove(Name_0_first)

    # Getting the second
    A_0_second = max(A_0)
    A_0_index = A_0.index(A_0_second)
    Name_0_second = Names_0[A_0_index]

    # Getting the first
    A_1_first = max(A_1)
    A_1_index = A_1.index(A_1_first)
    Name_1_first = Names_1[A_1_index]
    A_1.remove(A_1_first)
    Names_1.remove(Name_1_first)

    # Getting the second
    A_1_second = max(A_1)
    A_1_index = A_1.index(A_1_second)
    Name_1_second = Names_1[A_1_index]

    print("Strongest force at t = 0 and t = 0.5")
    print(f"|t = 0 : [{Name_0_first:10} : {A_0_first:7.4e} AU/Y^2]|, |t = 0.5 : [{Name_1_first:10} : {A_1_first:7.4e} AU/Y^2]|")
    print("Second strongest force at t = 0 and t = 0.5")
    print(f"|t = 0 : [{Name_0_second:10} : {A_0_second:7.4e} AU/Y^2]|, |t = 0.5 : [{Name_1_second:10} : {A_1_second:7.4e} AU/Y^2]|")

    print("--- Finished Gravity Tests ---\n\n")

    ax = plt.axes()
    ax.plot(r[0][0][int(T_0/dT):int(T_0/dT)+int(N*dt/dT)],r[1][0][int(T_0/dT):int(T_0/dT)+int(N*dt/dT)],label = "Planet 1 bane")
    ax.plot(ranged[0][1],ranged[1][1],label = "Planet 2 bane")
    ax.plot(R_T[0], R_T[1],color = "limegreen", label = "Simulert rakettbane")
        # Adding the star (not to scale)
    star = plt.Circle((0, 0), 0.75, color = 'gold')
    ax.add_patch(star)


    plt.xlabel("Position langs x-aksen [AU]")
    plt.ylabel("Position langs y-aksen [AU]")
    plt.title(f"Simulert rakettbane for t_0 = 2.2 Y, og t = {t} Y")
    plt.legend(loc = "upper right")
    plt.axis("equal")
    plt.show()


"""
Output:

---  Interpolation  Tests  ---
Start : [0 0], End : [1 0], Interpolator : 0.4, Expected Value : [0.4 0. ], Output : (0.4, 0.0)
Start : [2 0], End : [3 0], Interpolator : 0.4, Expected Value : [2.4 0. ], Output : (2.4, 0.0)
Start : [2 5], End : [3 4], Interpolator : 0.7, Expected Value : [2.7 4.3], Output : (2.7, 4.3)
Start : [2313.023 5283.382], End : [1000.2    1392.2321], Interpolator : 0.5, Expected Value : [1656.6115  3337.80705], Output : (1656.6115, 3337.80705)

Start : [0 0], End : [1 0], Interpolator : 0.4, Percent of diff : 0.40000
Start : [2 0], End : [3 0], Interpolator : 0.4, Percent of diff : 0.40000
Start : [2 5], End : [3 4], Interpolator : 0.7, Percent of diff : 0.70000
Start : [2313.023 5283.382], End : [1000.2    1392.2321], Interpolator : 0.5, Percent of diff : 0.50000
--- Finished Interpolation ---


---    Gravity    Tests    ---
Strongest force at t = 0 and t = 0.5
|t = 0 : [Planet 1   : 2.0273e+04 AU/Y^2]|, |t = 0.5 : [Sol        : 7.8854e+00 AU/Y^2]|
Second strongest force at t = 0 and t = 0.5
|t = 0 : [Sol        : 1.2167e+01 AU/Y^2]|, |t = 0.5 : [Planet 3   : 8.8837e-04 AU/Y^2]|
--- Finished Gravity Tests ---

"""