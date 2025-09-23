# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune


# Regular imports
import numpy             as np
import matplotlib.pyplot as plt

# Orbit imports
from SolarOrbit import SolarOrbit

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

if __name__ == "__main__":

    # This program creates a radial velocity curve with gaussian noise.
    # It starts by simulating the velocity, adding the peculiar velocity, which we decide, and then adding noise. 

    # Initializing AST
    Seed = utils.get_seed('bmthune')
    mission = SpaceMission(Seed)
    system = mission.system

    # Running the orbit simulation
    r_s,r_p,v_s,v_p,E = SolarOrbit(System= system,PlanetIndex = 2,dt = 1/10000, NumOrbits = 2)
    T = 2*np.pi*((system.semi_major_axes[2]**3)/(const.G_sol*(system.star_mass + system.masses[2])))**(1/2) * 2

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

    plt.title("Radiell hastighetskurve for stjernen")
    plt.xlabel("Tid [Y]")
    plt.ylabel("Radiell hastighet [AU/Y]")

    plt.show()
