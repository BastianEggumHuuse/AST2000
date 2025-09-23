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

def CM(r_1, m_1, r_2,m_2):

    m_2_2 = np.array([m_2, m_2]).T

    #CM_Raw = (m_1 * r_1 + m_2[0] *  r_2[0])/(m_1 + m_2[0])
    CM     = (m_1 * r_1 + np.sum(m_2_2  *  r_2,axis=0,keepdims=True))   /(m_1 + np.sum(m_2))

    return(CM)


def Orbit(r_1_0,r_2_0,v_1_0,v_2_0,m_1,m_2,T,dt):

    """
    Integration loop, performing a leapfrog integration on both bodies at once.

    r_1_0 : array(float) | Initial position for first body in AU (usually a star)
    r_2_0 : array(float) | Initial positions for all second bodies in AU (usually planets)
    v_1_0 : array(float) | Initial velocity for first body in AU/Y
    v_2_0 : array(float) | Initial velocitiess for all second bodies in AU/Y 
    m1    : float        | Mass of first body in solar masses
    m2    : float        | Masses of second bodies in solar masses
    T     : float        | Total time of simulation in Years
    dt    : float        | Timestep in years

    returns : 
    array(float) | Positional vectors for body 1
    array(float) | Positional vectors for all second bodies
    array(float) | Velocity vectors for body 1
    array(float) | Velocity vectors for all second bodies
    """

    # Initializing all arrays
    NumSteps = int(np.floor(T/dt))
    NumPlanets = len(r_2_0)
    r_1 = np.zeros((NumSteps,2))
    r_2 = np.zeros((NumSteps,NumPlanets,2))
    v_1 = np.zeros((NumSteps,2))
    v_2 = np.zeros((NumSteps,NumPlanets,2))
    a_1 = np.zeros((NumSteps,2))
    a_2 = np.zeros((NumSteps,NumPlanets,2))

    # We also initialize a center of mass array, which we use to adjust the cm every timestep
    # (technically this could just be a value, but we keep track for debugging purposes)
    cm = np.zeros((NumSteps,2))


    # Setting initial values for position and velocity
    r_1[0] = r_1_0
    r_2[0] = r_2_0
    v_1[0] = v_1_0
    v_2[0] = v_2_0

    # inner function that calculates the g-force from given planet at the given timestep
    def get_G(p,i):
        # Vector from body 1 to body 2
        r = r_2[i][p] - r_1[i]
        # Finding the unit of this vector
        r_hat = (r)/np.linalg.norm(r)
        # Calculating the G force, and multiplying with the unit to get a direction.
        G = (const.G_sol*m_1*m_2[p])/((np.linalg.norm(r))**2) *r_hat
        return G

    def multi_G(i):

        G = np.zeros((NumPlanets+1,2))
        for p in range(NumPlanets):
            G[p] = -(get_G(p,i))
        G[-1] = -(np.sum(G,axis = 0,keepdims = True))

        return(G)

    # Getting G force for the initial step
    G = multi_G(0)

    # Getting initial accelerations from G force
    a_1[0] = G[-1]/m_1
    a_2[0] = G[0:NumPlanets]/(np.array([m_2, m_2]).T)
    
    cm[0] = CM(r_1[0],m_1,r_2[0],m_2)

    # Running loop
    for i in range(0,NumSteps-1):

        # Leapfrog integrating the  (NOW VECTORIZED)
        # Note that we also subtract the CM! Since the CM moves, we have to adjust our frame every timestep.
        r_1[i+1] = r_1[i] + v_1[i]*dt + 0.5 * a_1[i]*dt**2
        r_2[i+1] = r_2[i] + v_2[i]*dt + 0.5 * a_2[i]*dt**2

        # Getting the g-force for the NEXT step
        G = multi_G(i+1)

        # j = 0
        # for g in G:
        #     j += 1
        #     print(f"{j}) {g}")
        # print("")
        # return

        # Getting the next acceleration
        a_1[i+1] = G[-1]/m_1
        a_2[i+1] = G[0:NumPlanets]/(np.array([m_2, m_2]).T)

        # Leapfrog integrating the velocity
        v_1[i+1] = v_1[i] + 0.5*(a_1[i] + a_1[i+1]) *dt
        v_2[i+1] = v_2[i] + 0.5*(a_2[i] + a_2[i+1]) *dt

        # Updating CM and Energy
        cm[i+1] = CM(r_1[i],m_1,r_2[i],m_2)
        
    # Returning all values
    return(r_1,r_2,v_1,v_2)
        
def SolarMultiOrbit(System,PlanetIndexes,dt,NumOrbits):

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
    m_p = System.masses[PlanetIndexes]
    # Setting initial positions, before changing frame
    r_s_0 = np.array([0,0])
    r_p_0 = System.initial_positions.T[PlanetIndexes]

    # Calculating center of mass
    cm = CM(r_s_0,m_s,r_p_0,m_p)

    # Moving our positions into the center of mass frame
    r_s = r_s_0 - cm
    r_p = r_p_0 - cm

    # Calculating total simulation time
    T = max(2*np.pi*((System.semi_major_axes[PlanetIndexes]**3)/(const.G_sol*(m_s + m_p)))**(1/2)) * NumOrbits

    # Initializing velocities, and running loop
    v_s_0 = np.array([0,0])
    v_p_0 = System.initial_velocities.T[PlanetIndexes] - cm
    r_s,r_p,v_s,v_p = Orbit(r_s,r_p,v_s_0,v_p_0,m_s,m_p,T,dt)

    # returning values
    return(r_s,r_p,v_s,v_p, T)

if __name__ == "__main__":

    # Initializing AST
    Seed = utils.get_seed('bmthune')
    mission = SpaceMission(Seed)
    system = mission.system
    
    r_s,r_p,v_s,v_p, T = SolarMultiOrbit(System= system,PlanetIndexes = [0,1,3,4,5],dt = 1/10000, NumOrbits = 1)
    
    # Plotting the two orbit
    # Initializing plotting
    fig, ax = plt.subplots()

    # Plotting orbits
    #ax.plot(r_s[:,0],r_s[:,1],color = "gold")
    # ax.plot(r_p.T[0][0],r_p.T[1][0],color = "violet")
    # ax.plot(r_p.T[0][1],r_p.T[1][1],color = "violet")
    # ax.plot(r_p.T[0][2],r_p.T[1][2],color = "violet")
    # ax.plot(r_p.T[0][3],r_p.T[1][3],color = "violet")

    # Velocity initialization
    v_pec = 0.1
    v_y = v_s[:,1] # Getting all y-axis velocities
    v_rad = v_y + v_pec 
    v_raw = v_rad.copy() # Getting a copy of the array pre-noise. When giving the array to another group, this will be commented out.

    # Gauss parameters
    mu = 0
    sigma = max(v_y)/5 

    # Adding noise to each timestep
    for i in range(len(v_rad)):
        v_rad[i] += np.random.normal(mu,sigma)
    
    # Plotting the noise graph and the no-noise graph
    plt.plot(np.linspace(0,T,len(v_rad)),v_rad)
    #plt.plot(np.linspace(0,T,len(v_rad)),v_raw,color = "red")

    v_noise = v_rad

    time = np.array([T])
    np.savez("SolarOrbitData",v_noise = v_noise,v_rad = v_rad,time = time)

    plt.title("Radiell hastighetskurve for vår stjerne")
    plt.xlabel("Tid [Y]")
    plt.ylabel("Radiell hastighet [AU/Y]")

    plt.show()