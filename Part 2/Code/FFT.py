import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import  ast2000tools.constants as const

seed = utils.get_seed('bmthune')
system = SolarSystem(seed)
"""""
print('My system has a {:g} solar mass star with a radius of {:g} kilometers.'
      .format(system.star_mass, system.star_radius))

for planet_idx in range(system.number_of_planets):
    print('Planet {:d} is a {} planet with a semi-major axis of {:g} AU and mass {:g}, radius = {:g}.'
          .format(planet_idx, system.types[planet_idx], system.semi_major_axes[planet_idx], system.masses[planet_idx]*1.989e+30/5.972e+24, system.radii[planet_idx]*1000/const.AU))


array = np.array([[1,1],[2,2],[3,3],[4,5],[6,7]])
array2 = np.array([[1,1],[2,2],[3,3],[4,2],[6,7]])
print(array/array2)
"""""
#system.print_info()




N = 2**10
T_tot = 10
sr = N/T_tot
n = np.linspace(0, T_tot, N)
x =  (
    1.0*np.sin(2*np.pi*1*n) +   # 1 Hz
    0.6*np.sin(2*np.pi*3*n) +   # 3 Hz
    0.3*np.sin(2*np.pi*7*n)     # 7 Hz
    
)
for i in range(len(x)):
      x[i] += np.random.normal(0,max(x)/5)

plt.plot(n,x)
plt.show()

def ditfft(X, N ):
      if N == 1:
            return X
            
      
      F = np.zeros(N, dtype=np.complex128)
      
      even = ditfft(X[:N:2], int(N/2))
      odd = ditfft(X[1:N:2], int(N/2))
      
      for k in range (0, N//2):
            p = even[k]
            q = np.exp(-1j *2 *np.pi/(N) *k)*odd[k]
            
            F[k] = (p + q)
            
            F[k + int(N/2)] = (p -q)
            
      return F
f = ditfft(x,N)
N = np.arange(len(f))
freq = N/T_tot
print(max(freq))
n_halv = len(N)//2
fre_halv = freq[:n_halv]
F_halv = f[:n_halv]


plt.plot(fre_halv,abs(F_halv))
plt.show()


    