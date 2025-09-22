# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune


# Regular imports
import numpy             as np
import matplotlib.pyplot as plt
from sys import exit

# Orbit imports
from OrbitPlotNumerical import NumericalOrbitFunction

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission



def Orbit(r_1_0,r_2_0,v_1_0,v_2_0,m_1,m_2,T,dt):

    NumSteps = int(np.floor(T/dt))
    r_1 = np.zeros((NumSteps,2))
    r_2 = np.zeros((NumSteps,2))
    v_1 = np.zeros((NumSteps,2))
    v_2 = np.zeros((NumSteps,2))
    a_1 = np.zeros((NumSteps,2))
    a_2 = np.zeros((NumSteps,2))
    E = np.zeros(NumSteps)

    cm = np.zeros((NumSteps,2))

    r_1[0] = r_1_0
    r_2[0] = r_2_0
    v_1[0] = v_1_0
    v_2[0] = v_2_0

    r = r_2[0] - r_1[0]
    r_hat = (r)/np.linalg.norm(r)
    G = (const.G_sol*m_1*m_2)/((np.linalg.norm(r))**2) *r_hat

    a_1[0] = G/m_1
    a_2[0] = G/m_2

    mu = (m_1*m_2)/(m_1 + m_2)
    E[0] = (0.5 * mu * np.linalg.norm(v_2[0] - v_1[0])**2) - (const.G_sol * m_1 * m_2)/np.linalg.norm(r_2[0] - r_1[0])

    t = 0
    for i in range(0,NumSteps-1):
        t =+ dt

        r_1[i+1] = r_1[i] + v_1[i]*dt + 0.5 * a_1[i]*dt**2 - cm[i]
        r_2[i+1] = r_2[i] + v_2[i]*dt + 0.5 * a_2[i]*dt**2 - cm[i]

        r = r_2[i+1] - r_1[i+1]
        r_hat = (r)/np.linalg.norm(r)
        G = (const.G_sol*m_1*m_2)/((np.linalg.norm(r))**2) *r_hat

        a_1[i+1] = G/m_1
        a_2[i+1] = -G/m_2

        v_1[i+1] = v_1[i] + 0.5*(a_1[i] + a_1[i+1]) *dt
        v_2[i+1] = v_2[i] + 0.5*(a_2[i] + a_2[i+1]) *dt

        cm[i+1] = (m_1 *  r_1[i] + m_2 *  r_2[i])/(m_1 + m_2)
        E[i+1] = (0.5 * mu * np.linalg.norm(v_2[i+1] - v_1[i+1])**2) - (const.G_sol * m_1 * m_2)/np.linalg.norm(r)

    return(r_1,r_2,v_1,v_2,E)
        


if __name__ == "__main__":

    # Initializing AST
    Seed = utils.get_seed('bmthune')
    mission = SpaceMission(Seed)
    system = mission.system

    # Planet 2 is a lot more massive than the other planets
    PlanetIndex = 2
    
    # Setting object masses
    m_s = system.star_mass
    m_p = system.masses[PlanetIndex]

    # Setting initial positions, before changing frame
    r_s_0 = np.array([0,0])
    r_p_0 = system.initial_positions.T[PlanetIndex]
    # Calculating center of mass
    cm = (m_s * r_s_0 + m_p * r_p_0)/(m_s + m_p)

    # Moving our positions into the center of mass frame
    r_s = r_s_0 - cm
    r_p = r_p_0 - cm

    # Calculating time
    dt = 1/10000
    NumOrbits = 2
    T = 2*np.pi*((system.semi_major_axes[PlanetIndex]**3)/(const.G_sol*(m_s + m_p)))**(1/2) * NumOrbits

    # Initializing velocities, and running loop
    v_s_0 = np.array([0,0])
    v_p_0 = system.initial_velocities.T[PlanetIndex] - cm
    r_s,r_p,v_s,v_p,E = Orbit(r_s,r_p,v_s_0,v_p_0,m_s,m_p,T,dt)

    # Output for C.1.2.b
    E_diff = (max(E)-min(E))/abs(min(E))
    print(f"Difference between maximum and minimum of E : {E_diff}")

    # Initializing plotting
    fig, ax = plt.subplots()

    # Plotting orbits
    ax.plot(r_s[:,0],r_s[:,1],color = "gold")
    ax.plot(r_p[:,0],r_p[:,1],color = "violet")

    # Adding title, and axis labels
    plt.title("Plott av planet og stjerne, i CM systemet")
    plt.xlabel("Posisjon langs x-aksen [AU]")
    plt.ylabel("Posisjon langs y-aksen [AU]")
    
    plt.show()
