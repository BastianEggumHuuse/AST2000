# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune


# Regular imports
import numpy             as np
import matplotlib.pyplot as plt
from sys import exit

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

def Orbit(r_1_0,r_2_0,v_1_0,v_2_0,m_1,m_2,T,dt):

    """
    Integration loop, performing a leapfrog integration on both bodies at once.

    r_1_0 : array(float) | Initial position for first body in AU (usually a star)
    r_2_0 : array(float) | Initial position for second body in AU (usually a planet)
    v_1_0 : array(float) | Initial velocity for first body in AU/Y
    v_2_0 : array(float) | Initial velocity for second body in AU/Y
    m1    : float        | Mass of first body in solar masses
    m2    : float        | Mass of second body in solar masses
    T     : float        | Total time of simulation in Years
    dt    : float        | Timestep in years

    returns : 
    array(float) | Positional vectors for body 1
    array(float) | Positional vectors for body 2
    array(float) | Velocity vectors for body 1
    array(float) | Velocity vectors for body 2
    array(float) | Energy at every timestep
    """

    # Initializing all arrays
    # They all have a dimention of (Numsteps,2)
    NumSteps = int(np.floor(T/dt))
    r_1 = np.zeros((NumSteps,2))
    r_2 = np.zeros((NumSteps,2))
    v_1 = np.zeros((NumSteps,2))
    v_2 = np.zeros((NumSteps,2))
    a_1 = np.zeros((NumSteps,2))
    a_2 = np.zeros((NumSteps,2))
    # Except for the energy of course :)
    E = np.zeros(NumSteps)

    # We also initialize a center of mass array, which we use to adjust the cm every timestep
    # (technically this could just be a value, but we keep track for debugging purposes)
    cm = np.zeros((NumSteps,2))

    # Setting initial values for position and velocity
    r_1[0] = r_1_0
    r_2[0] = r_2_0
    v_1[0] = v_1_0
    v_2[0] = v_2_0

    # inner function that calculates the g-force at the given timestep
    def get_G(i):
        # Vector from body 1 to body 2
        r = r_2[i] - r_1[i]
        # Finding the unit of this vector
        r_hat = (r)/np.linalg.norm(r)
        # Calculating the G force, and multiplying with the unit to get a direction.
        G = (const.G_sol*m_1*m_2)/((np.linalg.norm(r))**2) *r_hat
        return G

    # Getting G force for the initial step
    G = get_G(0)

    # Getting initial accelerations from G force
    a_1[0] = G/m_1
    a_2[0] = G/m_2

    # Defining the reduced mass
    mu = (m_1*m_2)/(m_1 + m_2)
    # Calculating the energy at the first step.
    E[0] = (0.5 * mu * np.linalg.norm(v_2[0] - v_1[0])**2) - (const.G_sol * m_1 * m_2)/np.linalg.norm(r_2[0] - r_1[0])

    # Running loop
    for i in range(0,NumSteps-1):

        # Leapfrog integrating the position
        # Note that we also subtract the CM! Since the CM moves, we have to adjust our frame every timestep.
        r_1[i+1] = r_1[i] + v_1[i]*dt + 0.5 * a_1[i]*dt**2 - cm[i]
        r_2[i+1] = r_2[i] + v_2[i]*dt + 0.5 * a_2[i]*dt**2 - cm[i]

        # Getting the g-force for the NEXT step
        G = get_G(i+1)
        # Getting the next acceleration
        a_1[i+1] = G/m_1
        a_2[i+1] = -G/m_2

        # Leapfrog integrating the velocity
        v_1[i+1] = v_1[i] + 0.5*(a_1[i] + a_1[i+1]) *dt
        v_2[i+1] = v_2[i] + 0.5*(a_2[i] + a_2[i+1]) *dt

        # Updating CM and Energy
        cm[i+1] = (m_1 *  r_1[i] + m_2 *  r_2[i])/(m_1 + m_2)
        E[i+1] = (0.5 * mu * np.linalg.norm(v_2[i+1] - v_1[i+1])**2) - (const.G_sol * m_1 * m_2)/np.linalg.norm(r_2[i+1] - r_1[i+1])

    # Returning all values
    return(r_1,r_2,v_1,v_2,E)
        
def SolarOrbit(System,PlanetIndex,dt,NumOrbits):

    """
    Wrapper for Orbit that takes care of all initialization

    PlanetIndex : int   | The index of the planet we want to simulate with the star
    dt          : float | Timestep in Years
    NumOrbits   : int   | Number of orbits we want to simulate, given in planet orbits (Note that these are calculated from a system where the star is static, so they aren't completely accurate)

    returns : 
    array(float) | Positional vectors for body 1
    array(float) | Positional vectors for body 2
    array(float) | Velocity vectors for body 1
    array(float) | Velocity vectors for body 2
    array(float) | Energy at every timestep
    """

    # Setting body masses
    m_s = System.star_mass
    m_p = System.masses[PlanetIndex]

    # Setting initial positions, before changing frame
    r_s_0 = np.array([0,0])
    r_p_0 = System.initial_positions.T[PlanetIndex]

    # Calculating center of mass
    cm = (m_s * r_s_0 + m_p * r_p_0)/(m_s + m_p)

    # Moving our positions into the center of mass frame
    r_s = r_s_0 - cm
    r_p = r_p_0 - cm

    # Calculating total simulation time
    T = 2*np.pi*((System.semi_major_axes[PlanetIndex]**3)/(const.G_sol*(m_s + m_p)))**(1/2) * NumOrbits

    # Initializing velocities, and running loop
    v_s_0 = np.array([0,0])
    v_p_0 = System.initial_velocities.T[PlanetIndex] - cm
    r_s,r_p,v_s,v_p,E = Orbit(r_s,r_p,v_s_0,v_p_0,m_s,m_p,T,dt)

    # returning values
    return(r_s,r_p,v_s,v_p,E)

if __name__ == "__main__":

    # Initializing AST
    Seed = utils.get_seed('bmthune')
    mission = SpaceMission(Seed)
    system = mission.system
    
    r_s,r_p,v_s,v_p,E = SolarOrbit(System= system,PlanetIndex = 2,dt = 1/10000, NumOrbits = 2)

    # Output for C.1.2.b
    E_diff = (max(E)-min(E))/abs(min(E))
    print(f"Difference between maximum and minimum of E : {E_diff}")

    # Plotting the two orvit
    # Initializing plotting
    fig, ax = plt.subplots()

    # Plotting orbits
    ax.plot(r_s[:,0],r_s[:,1],color = "gold")
    ax.plot(r_p[:,0],r_p[:,1],color = "violet")

    # Adding title, and axis labels
    plt.title("Plott av planet og stjerne, i CM systemet")
    plt.xlabel("Posisjon langs x-aksen [AU]")
    plt.ylabel("Posisjon langs y-aksen [AU]")
    
    # Making axes equal, and showing plot
    plt.axis('equal')
    plt.show()
