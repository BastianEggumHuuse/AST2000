import numpy             as np
import matplotlib.pyplot as plt
from numba import njit
import pickle as pkl 

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission
from HabitableZones import HabitableZones


#Unpickling launch
with open("Mission.pkl", 'rb') as file:#unplickeling the files
    mission = pkl.load(file)
#geting T_surface
Zones = HabitableZones(mission) 
Temp_s = Zones.Loop()[1]
T_s = np.mean(Temp_s)
#Mu gotten from previus task
mu = 44
#Setting other initial conditions and parameters for function
r_s = mission.system.radii[1]*1000
rho_s = mission.system.atmospheric_densities[1]
P_s = rho_s * const.k_B * T_s/(mu*const.m_p)
gamma = 1.4

#setting g to be the g halfway between ship and surface to get a better aproximation
g = const.G * mission.system.masses[1]*const.m_sun /(((r_s + 2684681)/2)**2)




def FindR(T, T_s):
    """
    A function to find what index we change zones by checking the temprature 
    Inputs T:Array of temprature values with same indexing as height and rho
           T_s: Surface temprature T_0
    Output: index of zoneswap 
    """
    for i in range(len(T)):
        if T[i] <= T_s/2:
            print(T[i])
            return i
    return len(T) 

def T_adi(rho):
    """
    Finding temprature in adiabatic zone, using formula given in the blogg
    Input Rho: Density, can be array
    Output: Temprature of same shape as Rho
    """
    a = P_s**(1-gamma) * T_s**gamma #finding the adiabatic constant bast on surface values
    T = a* (mu*const.m_p)**(1-gamma)/((const.k_B* rho)**(1-gamma))
    return T

def rho_abdi(r):
    """
    Finding the Rho in the adiabatic zone, using the formula given in the blogg. For certain R to high this will give an runtime warning du to taking root 
    of negative number, this is fine because this R will always be outside of the adibaatic zone and be overwritten later in the code
    inputs: r: can be array or number distance away from core not surface
    Outputs Rho in the same shape as R for all values of R
    """
    a = ( P_s**(1-gamma) * T_s**gamma)#finding the adiabatic constant bast on surface values
    C = 5/2 * rho_s**(2/5) +  g/(a*gamma) * (mu*const.m_p/const.k_B)**gamma * r_s
    rho = (2/5*(-g/(a*gamma) * (mu*const.m_p/const.k_B)**gamma * r + C))**(5/2)
    return rho

def TandRho(r):
    """
    Finding the total temprature and density over an array of r distnace from planet core
    Inputs r, distances from core (array)
    Outputs:
    T: Temprature array same shape as r
    rho: Density array same shape as r
    j: index for when zone changes 
    """
    rho = rho_abdi(r)#assumes adibatic for all r
    T = T_adi(rho)
    j = FindR(T, T_s) #finds when not adiabatic
    
   #overwrites all elements outside adiabatic area
    for i in range(j, len(T)):
        T[i] = T[j]
        rho[i] = rho[j]*np.exp(mu*const.m_p*g*(r[j]-r[i])/(const.k_B*T[i]))
    return T, rho, j
    




#setting an R array and finding density and temprature, array goes from surface to aprox rocket height        
R = np.linspace(r_s, 2.7e6, 100000)
T,Rho, j = TandRho(R)



#Doing a numerical integral over a column of air with base equal to 1m^2 going over every R   
Mass = 0
for i in range (len(Rho)-1):
   
    Mass += ((Rho[i] + Rho[i+1])/2 * (R[i+1]- R[i]))
    
#printing info dump
print(f'Masset av en søyle = {Mass}, Trykket ved overflaten over g = {P_s/g}')
print(f'Tethet ved overflate,    Model = {Rho[0] :.4} kg/m^3, Målt = {rho_s:.4} kg/m^3, relativ forskjel = {(Rho[0]-rho_s)/rho_s}')
print(f'Tempratur ved overflate, Model = {T[0] :.4} K       Målt = {T_s:.4} K     , relativ forskjel = {(T[0]-T_s)/T_s}')

#Plotting 
plt.vlines(R[j]-r_s, min(Rho),max(Rho), colors= 'red')
plt.plot(R-r_s,Rho)#minus R_s to get height over surface
plt.xlabel('Høyde [m]')
plt.ylabel('Tetthet [kg/m^3]')
plt.show()
plt.plot(R-r_s,T) #minus R_s to get height over surface
plt.vlines(R[j]-r_s, min(T),max(T), colors= 'red')
plt.xlabel('Høyde [m]')
plt.ylabel('Temp [K]')
plt.show()