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
with open("Mission.pkl", 'rb') as file:
    mission = pkl.load(file)

Zones = HabitableZones(mission)
Temp_s = Zones.Loop()[1]
T_s = np.mean(Temp_s)

mu = 16

r_s = mission.system.radii[1]*1000
rho_s = mission.system.atmospheric_densities[1]
P_s = rho_s * const.k_B * T_s/(mu*const.m_p)
gamma = 1.4


g = const.G * mission.system.masses[1]*const.m_sun /(r_s**2)


#for R i første del
def FindR(T, T_s):
    for i in range(len(T)):
        if T[i] <= T_s/2:
      
            return i
    return len(T) 

def T_adi(rho):
    
    a = P_s**(1-gamma) * T_s**gamma 

    T = a* (mu*const.m_p)**(1-gamma)/((const.k_B* rho)**(1-gamma))
    
    return T

def rho_abdi(r):
    a = P_s**(1-gamma) * T_s**gamma
    C = 5/2 * rho_s**(2/5) +  a*g/gamma * (mu*const.m_p/const.k_B)**gamma * r_s
    rho = (2/5*(-a*g/gamma * (mu*const.m_p/const.k_B)**gamma * r + C))**(5/2)
    return rho

def TandRho(r):
    rho = rho_abdi(r)
    T = T_adi(rho)
    j = FindR(T, T_s)
    
   
    for i in range(j, len(T)):
        T[i] = T[j]
        rho[i] = rho[j]*np.exp(mu*const.m_p*g/(const.k_B*T[i])*(r[j]-r[i]))
    return T, rho, j
    

def Rho(r,r_limit,T_limit):
    
    if r > r_limit:
        rho_limit = rho_abdi(r_limit)
        rho = rho_limit*np.exp(mu*const.m_p*g/(const.k_B*T_limit)*(r_limit-r))
    else:
        rho = rho_abdi(r)
    

    return rho

if __name__ == "__main__":

    R = np.linspace(r_s, 2305110.6806785692+1000, 1000000)

    T,Rhos, j = TandRho(R)
    #print(R[j], r_s)
    
    #print(T[j], T_s)

    Rhos = rho_abdi(R)
    plt.plot(R,Rhos)

    plt.show()
