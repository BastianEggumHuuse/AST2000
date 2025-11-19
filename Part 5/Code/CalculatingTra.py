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

# Importing the numerical planet positions and velocities
FilePath = "NumericalOrbitData.npz"
PlanetPositionFunction = NumericalOrbitFunction(FilePath)

# We want to use numba to speed up the computing process, but 
# NumericalOrbitFunction is a class, which numba isn't a huge fan of
# Therefore, we also load the data locally, to find planet positions.
npz = np.load(FilePath)
TotalTime    = float(npz["config"][0])
NumSteps     = int(npz["config"][2])
OrbitTimes   = npz["OrbitTimes"]
r            = npz["r"]
# Storing all this data in a handy tuple
Info = (TotalTime, NumSteps, OrbitTimes, r)

# Jit also doesn't like the constants class, so we save some constants we need.
AU = const.AU
G_sol = const.G_sol

@njit
def FindR(t,p, info):
    
    """
    Function that returns the position of a given planet along the x and y axes at a given time.
    This function is copied directly from the __call__ method of NumericalOrbitFunction,
    which is found within GeneralizedLaunch.py. See that file for comments.

    Parameters:
    t       : float        | the desired point in time
    p       : int          | the desired planet index
    info    : tuple        | tuple containing TotalTime,NumSteps,OrbitTimes, and planet positions.

    returns : Array(float) | the position of the given planet at the given time
    """
    
    TotalTime    = info[0]
    NumSteps     = info[1]
    OrbitTimes   = info[2]
    r            = info[3]
    
    if(t < 0):
        t = OrbitTimes[p] - t

    Index = int(np.floor((t/TotalTime)*NumSteps)) 

    x = (r[0][p][Index])
    y = (r[1][p][Index])
    return(np.array([x,y]))

@njit
def Lerp(R_0,R_1, I):

    """
    Linear interpolation function. Lerps between the 2D vectors R_0 and R_1.

    Parameters :
    R_0 : Array(float) | Vector we lerp from
    R_1 : Array(float) | Vector we lerp to
    I   : float        | How far between R_0 and R_1 to interpolate. Always between 0 and 1.

    returns :
    float | The x-coordinate of the interpolated vector
    float | The y-coordinate of the interpolated vector
    """

    # Finding the difference between R_0 and R_1
    dR_x = R_1[0] - R_0[0]
    dR_y = R_1[1] - R_0[1]
    
    # Lerping between R_0 and R_1 for both coordinates
    R_x = R_0[0] + I * dR_x
    R_y = R_0[1] + I * dR_y

    return R_x,R_y 

@njit
def GravitationalAks(R,K,N_k,dt, M,T_0, info):
    
    """
    Function that finds the total Gravitational Acceleration of the rocket.
    (We aquire a lot of parameters since numba doesn't like classes :[ )

    Parameters :
    R    : Array(float) | Array containing all positions
    K    : int          | The current timestep
    N_k  : int          | How many times the smaller timestep (rocket timestep) goes into the larger timestep (planet timestep)
    dt   : float        | Rocket timestep in years
    M    : Array(float) | Array containing all planet masses
    T_0  : float        | Start time in years (The time after reaching escape velocity)
    info : tuple        | Tuple containing info that FindR needs

    returns:
    float | total gravitational acceleration on the x-axis
    float | total gravitational acceleration on the y-axis
    """

    # Finding the gravitational acceleration from the sun (the suns mass is the last mass in the list)
    a_x = -G_sol *M[-1]*R[K][0]/ (((R[K][0])**2 + R[K][1]**2)**(3/2))
    a_y = -G_sol *M[-1]*R[K][1]/ (((R[K][0])**2 + R[K][1]**2)**(3/2))

    # Looping over all planets and computing gravitational acceleration
    for j in range(len(M)-1):
        
        # Finding positions that are closest to the current timestep
        # R_0 is backward in time, R_1 is forward in time
        R_0 = FindR(T_0 + K*dt, j, info)
        R_1 = FindR(T_0 + (K+N_k-1)*dt,j, info)

        # Computing the "true" planet position by lerping between R_0 and R_1
        R_p_x, R_p_y = Lerp(R_0, R_1, (K%N_k)/N_k) 

        # Computing difference between rocket position and planet position
        r_x = R_p_x - R[K][0] 
        r_y = R_p_y - R[K][1]

        # Computing most of the newtons gravitational equation
        # This bit is the same along both axes
        gamma = G_sol * M[j]/((r_x**2 + r_y**2)**(3/2))
            
        # Computing the gravitational acceleration on both axes, and adding them to the total acceleration
        a_x += r_x * gamma
        a_y += r_y * gamma

    return a_x, a_y

def GravitationAccelerations(R,K,N_k,dt, M,T_0, info):

    """
    This function does the exact same as the one above it, 
    but returns each acceleration individually (as opposed to the sum),
    as well as which bodies said accelerations come from.

    This function therefore has no comments :)
    """

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
def timestep(R,v,a, K, N_k, dt, M, T_0, info):

    """
    Function that computes one timestep of the trajectory, using leapfrog integration.

    Parameters :
    R    : Array(float) | Array containing all positions
    v    : Array(float) | Array containing all velocities
    a    : Array(float) | Array containing all accelerations
    K    : int          | The current timestep
    N_k  : int          | How many times the smaller timestep (rocket timestep) goes into the larger timestep (planet timestep)
    dt   : float        | Rocket timestep in years
    M    : Array(float) | Array containing all planet masses
    T_0  : float        | Start time in years (The time after reaching escape velocity)
    info : tuple        | Tuple containing info that FindR needs

    returns:
    None
    """

    R[K+1][0] = R[K][0] + v[K][0] * dt + 1/2*a[K][0]*dt**2
    R[K+1][1] = R[K][1] + v[K][1] * dt + 1/2*a[K][1]*dt**2

    a[K+1] = GravitationalAks(R,K+1,N_k,dt, M, T_0, info)
    
    v[K+1][0] = v[K][0] + 0.5*(a[K][0] + a[K+1][0])*dt
    v[K+1][1] = v[K][1] + 0.5*(a[K][1] + a[K+1][1])*dt

@njit
def Main(R_0, v_0, dt,dT, n, M,T_0, info):

    """
    The main function which computes the trajectory for all timesteps

    Parameters :
    R_0  : Array(float) | The initial position of the rocket
    v_0  : Array(float) | The initial velocity of the rocket
    dt   : float        | Rocket timestep in years
    dT   : float        | Planet timestep in years (larger than rocket)
    n    : int          | Amount of steps to simulate
    M    : Array(float) | Array containing all planet masses
    T_0  : float        | Start time in years (The time after reaching escape velocity)
    info : tuple        | Tuple containing info that FindR needs

    returns:
    Array(float) | Array containing all positions
    Array(float) | Array containing all velocities
    Array(float) | Array containing all accelerations
    float        | Time of simulation
    """

    # We simulate one more step than asked.
    # This is mainly for the case where we ask for 1 step, which normally would become 0 steps.
    N = n+1
    
    # Finding how many times dt goes into dT
    N_k = np.floor(dT/dt)
    
    # Defining our Arrays
    R = np.zeros((N,2))
    v = np.zeros((N,2))
    a = np.zeros((N,2))

    # Setting our initial values
    R[0] = R_0
    v[0] = v_0
    a[0] = np.array(GravitationalAks(R,0,N_k,dt, M,T_0, info))

    # Computing trajectory
    for K in range(N-1):
        timestep(R,v,a,K,N_k,dt,M,T_0, info)
        
    # Computing time
    T = N*dt

    return R,v,a, T

if __name__ == '__main__':
    
    # Interpolation Tests
    print("---  Interpolation  Tests  ---")

    # Testing multiple different interpolations against their expected outputs
    Interpolations = [
        [np.array([0,0]),np.array([1,0]),0.4,np.array([0.4,0])],
        [np.array([2,0]),np.array([3,0]),0.4,np.array([2.4,0])],
        [np.array([2,5]),np.array([3,4]),0.7,np.array([2.7,4.3])],
        [np.array([2313.023,5283.382]),np.array([1000.2,1392.2321]),0.5,np.array([1656.6115,3337.80705])],
        ]

    for i in Interpolations:
        print(f"Start : {i[0]}, End : {i[1]}, Interpolator : {i[2]}, Expected Value : {i[3]}, Output : {Lerp(i[0],i[1],i[2])}")

    print("")

    # Testing multiple different interpolations against their expected lengths
    for i in Interpolations:
        Lerped = np.array(Lerp(i[0],i[1],i[2])) - i[0]
        Lerped = np.linalg.norm(Lerped) / np.linalg.norm(i[1] - i[0])
        Lerped = np.linalg.norm(Lerped)
        print(f"Start : {i[0]}, End : {i[1]}, Interpolator : {i[2]}, Percent of diff : {Lerped:.5f}")

    print("--- Finished Interpolation --- \n\n")


    # Loading our spacemission class
    with open("Mission.pkl", 'rb') as file:
        mission = pkl.load(file)
    # Harvesting initial conditions
    R_0 = mission._position_after_launch
    v_0 =mission._velocity_after_launch
    M = np.zeros(mission.system._number_of_planets + 1)
    M[:-1] = mission.system.masses
    M[-1] = mission.system.star_mass
    T_0 = mission.time_after_launch
    dT = PlanetPositionFunction.dt
    dt = 1/1000000
    T_0 = mission.time_after_launch

    # Defining length of simulation, and translating to rocket timesteps
    t = 3.14
    N = int(t / dt)#3140000 * 2
    
    # Running simulation
    R, v, a, T = Main(R_0, v_0, dt,dT, N, M ,T_0, Info)
    # Transposing array to access our data
    R_T = R.T

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

    # Printing the two strongest forces at t = 0Y and t = 0.5Y (after launch)
    print("Strongest force at t = 0 and t = 0.5")
    print(f"|t = 0 : [{Name_0_first:10} : {A_0_first:7.4e} AU/Y^2]|, |t = 0.5 : [{Name_1_first:10} : {A_1_first:7.4e} AU/Y^2]|")
    print("Second strongest force at t = 0 and t = 0.5")
    print(f"|t = 0 : [{Name_0_second:10} : {A_0_second:7.4e} AU/Y^2]|, |t = 0.5 : [{Name_1_second:10} : {A_1_second:7.4e} AU/Y^2]|")

    print("--- Finished Gravity Tests ---\n\n")

    # Plotting the trajectory along with the trajectory of the planets.
    ax = plt.axes()
    ax.plot(r[0][0][int(T_0/dT):int(T_0/dT)+int(N*dt/dT)],r[1][0][int(T_0/dT):int(T_0/dT)+int(N*dt/dT)],label = "Planet 1 bane")
    ranged = PlanetPositionFunction.range(0,6)
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