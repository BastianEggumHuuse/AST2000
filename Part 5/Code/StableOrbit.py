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

# Constants
G = const.G

def OrbitSimulation(spacecraft_pos, spacecraft_vel, spacecraft_mass, planet_pos,planet_vel, planet_mass):

    """
    This function performs the entire orbit simulation, displays the results, and plots the final simulated orbit.

    Parameters:
    spacecraft_pos  : Array(float) | Position of spacecraft in solar system frame (in m)
    spacecraft_vel  : Array(float) | Velocity of spacecraft in solar system frame (in m/s)
    spacecraft_mass : float        | Mass of spacecraft (in kg)
    planet_pos      : Array(float) | Position of planet in solar system frame (in m)
    planet_vel      : Array(float) | Velocity of planet in solar system frame (in m/s)
    planet_mass     : float        | Mass of planet (in kg)
    """

    """
    Part 1, finding the initial conditions
    """

    # Finding distance from planet
    r_vec = spacecraft_pos - planet_pos
    r = np.linalg.norm(r_vec)

    # Finding radial velocity
    r_hat = r_vec / np.linalg.norm(r_vec)
    relative_vel = spacecraft_vel - planet_vel
    v_r = np.dot(relative_vel,r_hat)

    # Finding angle
    if(r_vec[0] == 0):
        # Not included in the flowchart.
        # If x-coordinate is 0, this division is illegal, so we introduce a small number instead
        theta = np.arctan(r_vec[1]/(0.000001))
    else:
        theta = np.arctan(r_vec[1]/r_vec[0])
    if(r_vec[0] < 0):
        theta += (r_vec[1]/abs(r_vec[1]))*np.pi

    # Finding angular velocity
    v_theta_vec = relative_vel - (v_r * r_hat) # v_theta is orthogonal to v_r, so we simply remove the part of Spacecraft_vel that corresponds to the component v_r
    v_theta = np.linalg.norm(v_theta_vec)/r

    # Finding the direction of the angular velocity
    dir = np.cross(r_hat,v_theta_vec)/np.linalg.norm(np.cross(v_theta_vec,r_hat))
    v_theta *= dir

    # Finding the angular momentum (using v_theta and r)
    # The angular momentum is conserved through the entire movement.
    l = ReducedMass(spacecraft_mass,planet_mass)*(r**2)*v_theta

    """
    Part 2, simulating the two-body-problem
    """

    # Performing a singular orbit
    R,V_R,THETA,V_THETA,t_0 = SingleOrbit(r,v_r,theta,v_theta,spacecraft_mass,planet_mass,l)

    # Performing analysis
    print("Printing data from first analysis:\n")
    peri_0,apo_0,a_0,b_0,e_0 = OrbitAnalysis(R,THETA)
    PrintAnalysis(a_0,b_0,e_0,t_0,apo_0,peri_0)

    # Performing 2 more orbits
    for i in range(2):
        R,V_R,THETA,V_THETA,t = SingleOrbit(R[-1],V_R[-1],THETA[-1],V_THETA[-1],spacecraft_mass,planet_mass,l)

    # Performing some more analysis
    print("Printing data from second analysis:\n")
    peri,apo,a,b,e = OrbitAnalysis(R,THETA)
    PrintAnalysis(a,b,e,t,apo,peri)

    print("Printing relative differences between first and second analyses:\n")
    DiffAnalysis((a_0,a),(b_0,b),(e_0,e),(t_0,t),(apo_0,apo),(peri_0,peri))

    # Finding out how many orbits are in two years, using the previously calculated orbit time
    two_years = 2*(60*60*24*365)
    num_orbits = int(np.floor(two_years/t))

    # Performing two year's worth of orbits
    total_time = 0
    for i in range(num_orbits):
        R,V_R,THETA,V_THETA,t = SingleOrbit(R[-1],V_R[-1],THETA[-1],V_THETA[-1],spacecraft_mass,planet_mass,l)
        total_time += t

    # Performing the final analysis
    print(f"Printing data from final analysis conducted after {(total_time * 2) /(two_years)} years.:\n")
    peri,apo,a,b,e = OrbitAnalysis(R,THETA)
    PrintAnalysis(a,b,e,t,apo,peri)

    print("Printing relative differences between first and final analyses:\n")
    DiffAnalysis((a_0,a),(b_0,b),(e_0,e),(t_0,t),(apo_0,apo),(peri_0,peri))

    # Plotting final orbit
    plt.plot((planet_pos[0] + R*np.cos(THETA)),(planet_pos[1] + R*np.sin(THETA)))
    plt.plot(planet_pos[0],planet_pos[1],"o")
    plt.axis("equal")

    plt.show()

@njit
def SingleOrbit(r,v_r,theta,v_theta,m1,m2,l):

    """
    Function that performs a single orbit around the planet

    Parameters:
    r       : float | initial distance between planet and spacecraft
    v_r     : float | initial radial velocity between planet and spacecraft
    theta   : float | initial angle between planet and spacecraft
    v_theta : float | initial angular velocity between planet and spacecraft
    m1      : float | mass of spacecraft
    m2      : float | mass of planet
    l       : float | angular momentum for the system (this is a conserved value)

    returns:
    Array(float) | All calculated distances between planet and spacecraft
    Array(float) | All calculated radial velocities between planet and spacecraft
    Array(float) | All calculated angles between planet and spacecraft
    Array(float) | All calculated angular velocities between planet and spacecraft
    float        | total orbit time.
    """

    # Defining our time interval, and our maximum time
    # We assume that the orbit will take no longer that 10 days to complete
    dt = 1
    t = 10*60*60*24
    N = int(t/dt)

    # Defining our arrays
    R = np.zeros(N)
    V_R = np.zeros(N)
    THETA = np.zeros(N)
    V_THETA = np.zeros(N)

    # Setting initial values
    R[0] = r
    V_R[0] = v_r
    THETA[0] = theta
    V_THETA[0] = v_theta

    # Defining acceleration array (this is not returned, but neccesary for leapfrog integration)
    A_R = np.zeros(N)
    A_R[0] = RadialAcceleration(R[0],V_THETA[0],m1,m2)

    # Looping
    for i in range(1,N):

        # Updating R
        R[i] = R[i-1] + V_R[i-1] * dt + 0.5 * A_R[i-1] * dt**2
        # Updating A_R
        A_R[i] =  RadialAcceleration(R[i],V_THETA[i-1],m1,m2)
        # Updating V_R 
        V_R[i] = V_R[i-1] + 0.5*(A_R[i-1] + A_R[i]) * dt

        # Updating V_THETA
        V_THETA[i] = AngularVelocity(R[i],m1,m2,l) * dt
        # Updating THETA
        THETA[i] = THETA[i-1] + V_THETA[i] * dt

        # Checking if we have completed a loop
        if(abs(THETA[i]) > 2*np.pi + abs(theta)):
            break

    return R[:i],V_R[:i],THETA[:i],V_THETA[:i],i*dt

def OrbitAnalysis(R,THETA):

    """
    Function that analyses the data from an orbit and computes several attributes of the orbit

    Parameters : 
    R     : Array(float) | All computed distances between planet and spacecraft for one orbit
    THETA : Array(float) | All computed angles between planet and spacecraft for one orbit

    Returns :
    float | periapsis of the orbit
    float | apoapsis of the orbit
    float | semi-major axis of the orbit
    float | semi-minor axis of the orbit
    float | eccentricity of the orbit
    """

    periapsis = min(R)
    apoapsis = max(R)
    semi_major = (periapsis + apoapsis)/2
    eccentricity = (semi_major - periapsis)/semi_major
    semi_minor = ((semi_major**2)*(1-eccentricity**2))**0.5

    return periapsis,apoapsis,semi_major,semi_minor,eccentricity

def PrintAnalysis(a,b,e,p,apo,peri):

    """
    Function that prints the data returned by OrbitAnalysis (along with time)
    Note that distances are converted from meters to kilometers,
    and time intervals are converted from seconds to hours
    """

    print(f"Semi-major axis : {a / 1000:10.3f} km")
    print(f"Semi-minor axis : {b / 1000:10.3f} km")
    print(f"Eccentricity    : {e:10.3f}")
    print(f"Orbital period  : {p/60/60:10.3f} hours")
    print(f"Apoapsis        : {apo/1000:10.3f} km")
    print(f"Periapsis       : {peri/1000:10.3f} km\n")

def DiffAnalysis(a,b,e,p,apo,peri):

    """
    Function that prints the relative differences between two sets of data returned by OrbitAnalysis (along with time)
    All parameters are tuples containing two values from to different orbits.
    """

    print(f"Semi-major axis : {RelativeTuple(a):.7f}")
    print(f"Semi-minor axis : {RelativeTuple(b):.7f}")
    print(f"Eccentricity    : {RelativeTuple(e):.7f}")
    print(f"Orbital period  : {RelativeTuple(p):.7f}")
    print(f"Apoapsis        : {RelativeTuple(apo):.7f}")
    print(f"Periapsis       : {RelativeTuple(peri):.7f}\n")

def RelativeTuple(t):
    # Function that computes the relative difference between two values contained in a tuple
    return abs((t[0]-t[1])/t[1])

@njit
def GravitationalForce(r,m1,m2):
    # Function that computes the gravitational force between two objects
    return (G*m1*m2/r**2)

@njit
def ReducedMass(m1,m2):
    # Function that computes the reduced mass between two objects
    return (m1*m2)/(m1+m2)

@njit
def RadialAcceleration(r,v_theta,m1,m2):
    # Function that computes the radial acceleration between two objects in a two-body system.
    a_r = -GravitationalForce(r,m1,m2)/ReducedMass(m1,m2) + r*v_theta**2
    return a_r

@njit
def AngularVelocity(r,m1,m2,l):
    # Function that computes the angular velocity between two objects in a two-body system.
    v_theta = l/(ReducedMass(m1,m2)*r**2)
    return v_theta


if __name__ == "__main__":

    # Ast init
    seed = utils.get_seed('bmthune')
    mission = SpaceMission(seed)

    Filepath = "NumericalOrbitData.npz"
    planet_positions = NumericalOrbitFunction(Filepath)

    r_0 = np.array((2.91050526,2.50762015)) * const.AU
    v_0 = np.array((-3.44042026,3.54089592)) * (const.AU / (60*60*24*365))
    t_0 = 3.8499120193480745
    m_0 = mission.spacecraft_mass + 37.4

    r_p = planet_positions(t_0,1) * const.AU
    v_p = planet_positions.GetVelocity(t_0,1) * (const.AU / (60*60*24*365))
    m_p = mission.system.masses[1] * const.m_sun

    #v_0 = v_0 - v_p

    # Printing initial conditions
    print("Printing inital conditions: ")
    print(f"Time of simulation: {t_0:7.4e}")
    print(f"Rocket : |r_0 : [{r_0[0]:7.5e},{r_0[1]:7.5e}], v_0 : [{v_0[0]:7.5e},{v_0[1]:7.5e}], m_0 : {m_0:7.5e}|")
    print(f"Planet : |r_p : [{r_p[0]:7.5e},{r_p[1]:7.5e}], v_p : [{v_p[0]:7.5e},{v_p[1]:7.5e}], m_p : {m_p:7.5e}|\n\n")

    # Displaying initial conditions
    r_l = r_0 - r_p
    v_l = v_0 - v_p
    ax = plt.axes()
    # Adding the rocket velocity
    ax.quiver(r_l[0],r_l[1],v_l[0],v_l[1])
    # Adding the rocket
    ax.scatter(r_l[0],r_l[1],color = "firebrick")
    # Adding the planet
    planet = plt.Circle((0, 0), r_l[0]/10, color = 'royalblue')
    ax.add_patch(planet)
    # Showing plot
    plt.xlabel("Posisjon langs x-aksen [m]") 
    plt.ylabel("Posisjon langs y-aksen [m]") 
    plt.axis("equal")
    plt.show()

    OrbitSimulation(r_0,v_0,m_0,r_p,v_p,m_p)

