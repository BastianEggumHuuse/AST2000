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

# Assume a Spacecraft position, a Planet position, and a Spacecraft velocity
def OrbitSimulation(spacecraft_pos, spacecraft_vel, spacecraft_mass, planet_pos, planet_mass):

    """
    Part 1, finding the initial conditions
    """

    # Finding distance from planet
    r_vec = spacecraft_pos - planet_pos
    r = np.linalg.norm(r_vec)

    # Finding radial velocity
    r_hat = r_vec / np.linalg.norm(r_vec)
    v_r = np.dot(spacecraft_vel,r_hat)

    # Finding angle
    if(r_vec[0] == 0):
        theta = np.arctan(r_vec[1]/(0.000001))
    else:
        theta = np.arctan(r_vec[1]/r_vec[0])
    if(r_vec[0] < 0):
        theta += (r_vec[1]/abs(r_vec[1]))*np.pi

    # Finding angular velocity
    v_theta_vec = spacecraft_vel - (v_r *r_hat) # v_theta is orthogonal to v_r, so we simply remove the part of Spacecraft_vel that corresponds to the component v_r
    v_theta = np.linalg.norm(v_theta_vec)/r

    # Finding the direction of the angular velocity
    dir = np.cross(r_hat,v_theta_vec)/np.linalg.norm(np.cross(v_theta_vec,r_hat))
    v_theta *= dir

    # Finding the angular momentum (using v_theta and r)
    l = ReducedMass(spacecraft_mass,planet_mass)*(r**2)*v_theta

    """
    Part 2, simulating the two-body-problem
    """

    # Performing a singular orbit
    R,V_R,THETA,V_THETA,t = SingleOrbit(r,v_r,theta,v_theta,spacecraft_mass,planet_mass,l)

    # Performing analysis
    print("Printing data from first analysis:\n")
    peri,apo,a,b,e = OrbitAnalysis(R,THETA)
    PrintAnalysis(a,b,e,t,apo,peri)

    # Performing 2 more loops
    for i in range(2):
        R,V_R,THETA,V_THETA,t = SingleOrbit(R[-1],V_R[-1],THETA[-1],V_THETA[-1],spacecraft_mass,planet_mass,l)

    # Performing some more analysis
    print("Printing data from second analysis:\n")
    peri,apo,a,b,e = OrbitAnalysis(R,THETA)
    PrintAnalysis(a,b,e,t,apo,peri)

    # Plotting final orbit
    plt.plot((planet_pos[0] + R*np.cos(THETA)),(planet_pos[1] + R*np.sin(THETA)))
    plt.plot(planet_pos[0],planet_pos[1],"o")
    plt.axis("equal")

    plt.show()

@njit
def SingleOrbit(r,v_r,theta,v_theta,m1,m2,l):

    dt = 1
    t = 60*60*24
    N = int(t/dt)

    R = np.zeros(N)
    V_R = np.zeros(N)
    THETA = np.zeros(N)
    V_THETA = np.zeros(N)

    R[0] = r
    V_R[0] = v_r
    THETA[0] = theta
    V_THETA[0] = v_theta

    A_R = np.zeros(N)
    A_R[0] = RadialAcceleration(R[0],V_THETA[0],m1,m2)

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

    # Finding periapsis, apoapsis, and semi_major_axis
    periapsis = min(R)
    apoapsis = max(R)
    semi_major = periapsis + apoapsis

    # Finding semi minor requires a loop:
    peri_index = int(np.where(R == periapsis)[0][0])
    peri_angle = THETA[peri_index]
    semi_minor_angle = peri_angle + np.pi/2
    #semi_minor_index = np.where()

    for i in range(peri_index,len(THETA)-1):
        if(THETA[i] % semi_minor_angle < THETA[i+1] % semi_minor_angle):
            break

    semi_minor = R[i] * 2
    eccentricity = (1-(semi_minor/semi_major)**2)**0.5

    return periapsis,apoapsis,semi_major,semi_minor,eccentricity

def PrintAnalysis(a,b,e,p,apo,peri):

    print(f"Semi-major axis : {a:20} m")
    print(f"Semi-minor axis : {b:20} m")
    print(f"Eccentricity    : {e:20}")
    print(f"Orbital period  : {p:20} s")
    print(f"Apoapsis        : {apo:20} m")
    print(f"Apoapsis        : {peri:20} m\n")

@njit
def GravitationalForce(r,m1,m2):
    return (G*m1*m2/r**2)

@njit
def ReducedMass(m1,m2):
    return (m1*m2)/(m1+m2)

@njit
def RadialAcceleration(r,v_theta,m1,m2):
    a_r = -GravitationalForce(r,m1,m2)/ReducedMass(m1,m2) + r*v_theta**2
    return a_r

@njit
def AngularVelocity(r,m1,m2,l):
    v_theta = l/(ReducedMass(m1,m2)*r**2)
    return v_theta

if __name__ == "__main__":

    # Ast init
    seed = utils.get_seed('bmthune')
    mission = SpaceMission(seed)

    Filepath = "NumericalOrbitData.npz"
    planet_positions = NumericalOrbitFunction(Filepath)

    t_0 = 3

    # Planet info
    planet_position = planet_positions(3,1) * const.AU
    planet_radius = mission.system.radii[1] * 1000
    planet_mass = mission.system.masses[1] * const.m_sun

    # Spacecraft info
    spacecraft_direction = np.array([1,0])/np.linalg.norm(np.array([1,0]))
    spacecraft_position = planet_position + spacecraft_direction * planet_radius * 3
    spacecraft_velocity = np.array([0.1,1]) * ((G*planet_mass)/np.linalg.norm(spacecraft_position - planet_position))**0.5
    spacecraft_mass = mission.spacecraft_mass

    OrbitSimulation(spacecraft_position,spacecraft_velocity,spacecraft_mass,planet_position,planet_mass)

