# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune


# Regular imports
import numpy             as np
import matplotlib.pyplot as plt

# Orbit imports
from SolarOrbit import Orbit

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission


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
    NumOrbits = 1
    T = 2*np.pi*((system.semi_major_axes[PlanetIndex]**3)/(const.G_sol*(m_s + m_p)))**(1/2) * NumOrbits

    # Initializing velocities, and running loop
    v_s_0 = np.array([0,0])
    v_p_0 = system.initial_velocities.T[PlanetIndex] - cm
    r_s,r_p,v_s,v_p,E = Orbit(r_s,r_p,v_s_0,v_p_0,m_s,m_p,T,dt)

    
    # Game start
    v_pec = 0.1 # men bendik er crazy
    v_y = v_s[:,1]
    v_rad = v_y + v_pec
    v_backup = v_rad.copy()
    mu = 0
    sigma = max(v_y)/5 

    for i in range(len(v_rad)):
        v_rad[i] += np.random.normal(mu,sigma)
    
    plt.plot(np.linspace(0,T,len(v_rad)),v_rad)
    plt.plot(np.linspace(0,T,len(v_rad)),v_backup,color = "red")

    plt.title("Radiell hastighetskurve for vår stjerne")
    plt.xlabel("Tid [AU]")
    plt.ylabel("Hastighet [AU/Y]")

    plt.show()
