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

def Area_Circle(r): #defining fast way to take the area of a circle, to save time
    return(np.pi * r **2)

if __name__ == "__main__":

    # Initializing AST 
    Seed = utils.get_seed('bmthune')
    mission = SpaceMission(Seed)
    system = mission.system

    # Running the orbit simulation
    dt = 1/100000 #setting delta time equal to the one used in previus codes
    r_s,r_p,v_s,v_p,E = SolarOrbit(System= system,PlanetIndex = 2,dt = dt, NumOrbits = 1) #Running solar orbits to get the vellocity of planet and star
    T = 2*np.pi*((system.semi_major_axes[2]**3)/(const.G_sol*(system.star_mass + system.masses[2])))**(1/2) #Calculating the orbit time of one rotation from keeplers third law
    Radius_s = (system.star_radius * 1000) / const.AU;Radius_p = (system.radii[2] * 1000) / const.AU #Fidning the radius of the sun and the planet in AU

    F = np.ones(len(r_s)) #Creating the relative flux array, should be one whenever no planet in the way
    for i in range(len(r_s)):   #Looping over all the timesteps for one rotation. 

        r_s_x = r_s[i][0] #fidning the x cordinate of both the stare and planet
        r_p_x = r_p[i][0]

        # Case 1 Checking if the planet is in front of the star, cheking if the radi are crossing eachother
        if not ((r_s_x - Radius_s - Radius_p < r_p_x) and (r_p_x < r_s_x + Radius_s + Radius_p ) and r_p[i][1] < r_s[i][1]):
            continue
        
        #creating the vector from sun to plannet
        r_sp =  (r_p_x - r_s_x) 
        #creating the normalized vector of the one over
        r_sp_hat = (r_sp/abs(r_sp))
        #Defining r_o from the planets core to the sun radi, multiplying it with the dir vector to make sure we take both the way inn and out into acount
        r_o = r_sp_hat*(r_sp - (r_s_x + (Radius_s*r_sp_hat)))

        # Case 2 When we are enterly on the inside of the star, just sett the relative flux to the 1 - Ap/As- 
        if (r_s_x - Radius_s + Radius_p < r_p_x) and (r_p_x < r_s_x + Radius_s - Radius_p ):
            F[i] = (Area_Circle(Radius_s) - Area_Circle(Radius_p))/Area_Circle(Radius_s)

            continue

        # Case 3 this is when its on its way in or out
        #finding the angle between the x axies the planet center and where the star crosses the planet, assumes the star is a straight line
        theta = np.arccos(r_o/Radius_p)
        
        #finding the slice of the sircle between the slicing points
        dA = (Radius_p**2 * theta)
        #Fidning the area of the part inside the star, doing this by subtracting the triangle on the otherside of the line from the pie (se float chart)
        dA_s = dA - Radius_p * np.sin(theta) * r_o

        #Seting the relative flux equal to the the 1 - the relation between the area crossing the sun and the total are of the sun 
        F[i] = (Area_Circle(Radius_s) - dA_s)/Area_Circle(Radius_s)
        
        #Adding some gausian noise
    for i in range(len(F)):
         F[i] += np.random.normal(0,1e-4)

    # Plotting the noise graph
    #making arrays for plotting
    t_0 = int(np.floor(18.882/dt))
    t_1 = int(np.floor(18.892/dt))
    plt.plot(np.linspace(t_0,t_1,len(F[t_0:t_1])),F[t_0:t_1])
    #adding title and lables
    plt.title("Lyskurve for vår stjerne")
    plt.xlabel("Tid [Y]")
    plt.ylabel("Radiell hastighet [AU/Y]")

    plt.show()

    

            