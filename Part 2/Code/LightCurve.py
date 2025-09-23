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

def Area_Circle(r):
    return(np.pi * r **2)

if __name__ == "__main__":

    # Initializing AST
    Seed = utils.get_seed('bmthune')
    mission = SpaceMission(Seed)
    system = mission.system

    # Running the orbit simulation
    r_s,r_p,v_s,v_p,E = SolarOrbit(System= system,PlanetIndex = 2,dt = 1/10000, NumOrbits = 2)
    T = 2*np.pi*((system.semi_major_axes[2]**3)/(const.G_sol*(system.star_mass + system.masses[2])))**(1/2) * 2
    Radius_s = (system.star_radius * 1000) / 149597870700;Radius_p = (system.radii[2] * 1000) / 149597870700

    F = np.zeros(len(r_s))
    for i in range(len(r_s)):

        r_s_x = r_s[i][0]
        r_p_x = r_p[i][0]

        # Case 1
        if not ((r_s_x - Radius_s - Radius_p < r_p_x) and (r_p_x < r_s_x + Radius_s + Radius_p)):
            F[i] = 1
            continue

        r_o = r_p_x - 2*r_s_x

        # Case 2
        if(abs(r_o) > Radius_p):
            F[i] = (Area_Circle(Radius_s) - Area_Circle(Radius_p))/Area_Circle(Radius_s)
            if(F[i] > 1):
                print(f"Inner : {F[i]}")
            continue

        # Case 3
        theta = 2 * np.arccos(r_o/Radius_p)
        
        dA = (Radius_p**2 * theta)/2
        dA_s = dA - Radius_p * np.sin(theta) * r_o

        F[i] = (Area_Circle(Radius_s) - dA_s)/Area_Circle(Radius_s)
        if(F[i] > 1):
                print(f"Transition : {dA,Radius_p * np.sin(theta) * r_o}")

    # Plotting the noise graph and the no-noise graph
    plt.plot(np.linspace(0,T,len(F)),F)

    plt.title("Radiell hastighetskurve for vår stjerne")
    plt.xlabel("Tid [AU]")
    plt.ylabel("Radiell hastighet [AU/Y]")

    plt.show()

    

            