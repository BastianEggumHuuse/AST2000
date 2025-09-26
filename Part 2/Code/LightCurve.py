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
    dt = 1/100000
    r_s,r_p,v_s,v_p,E = SolarOrbit(System= system,PlanetIndex = 2,dt = dt, NumOrbits = 1)
    T = 2*np.pi*((system.semi_major_axes[2]**3)/(const.G_sol*(system.star_mass + system.masses[2])))**(1/2) * 1
    Radius_s = (system.star_radius * 1000) / const.AU;Radius_p = (system.radii[2] * 1000) / const.AU

    F = np.ones(len(r_s))
    for i in range(len(r_s)):   

        r_s_x = r_s[i][0]
        r_p_x = r_p[i][0]

        # Case 1
        if not ((r_s_x - Radius_s - Radius_p < r_p_x) and (r_p_x < r_s_x + Radius_s + Radius_p ) and r_p[i][1] < r_s[i][1]):
            continue

        r_sp =  (r_p_x - r_s_x) 
        r_sp_hat = (r_sp/abs(r_sp))
        r_o = r_sp_hat*(r_sp - (r_s_x + (Radius_s*r_sp_hat)))

        # Case 2
        if (r_s_x - Radius_s + Radius_p < r_p_x) and (r_p_x < r_s_x + Radius_s - Radius_p ):
            F[i] = (Area_Circle(Radius_s) - Area_Circle(Radius_p))/Area_Circle(Radius_s)
            if(F[i] > 1):
                print(f"Inner : {F[i]}")
            continue

        # Case 3
        theta = np.arccos(r_o/Radius_p)
        
        dA = (Radius_p**2 * theta)
        dA_s = dA - Radius_p * np.sin(theta) * r_o

        F[i] = (Area_Circle(Radius_s) - dA_s)/Area_Circle(Radius_s)
        if(F[i] > 1):
                print(f"Transition : {dA,Radius_p * np.sin(theta) * r_o}")

    for i in range(len(F)):
         F[i] += np.random.normal(0,1e-4)

    # Plotting the noise graph
    t_0 = int(np.floor(18.882/dt))
    t_1 = int(np.floor(18.892/dt))
    plt.plot(np.linspace(t_0,t_1,len(F[t_0:t_1])),F[t_0:t_1])

    plt.title("Lyskurve for vår stjerne")
    plt.xlabel("Tid [Y]")
    plt.ylabel("Radiell hastighet [AU/Y]")

    plt.show()

    

            