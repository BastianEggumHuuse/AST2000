import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import  ast2000tools.constants as const

seed = utils.get_seed('bmthune')
system = SolarSystem(seed)

def ditfft(X, N):
      if N == 1:
            return X


      F = np.zeros(N, dtype=np.complex128)

      even = ditfft(X[:N:2], int(N/2))
      odd = ditfft(X[1:N:2], int(N/2))

      for k in range (0, N//2):
            p = even[k]
            q = np.exp(-1j * 2*np.pi/(N) * k) * odd[k]

            F[k] = (p + q)

            F[k + int(N/2)] = (p -q)

      return F

FileName = "SolarOrbitData.npz"
RawData = np.load(FileName)

TotalTime = RawData['time'][0]
N = 2**20

Data = RawData['v_raw'][:int(np.floor(430/TotalTime *len(RawData['v_raw'])))]
TotalTime = 430
v_pek = np.mean(Data)



NoiseData = Data -v_pek

TimeSample = np.linspace(0,TotalTime,N)
dt = TimeSample[1]- TimeSample[0]

index = np.linspace(0, len(NoiseData)-1,N, dtype=int)
NoiseSample = NoiseData[index] 



f = ditfft(NoiseSample,N)
f_range = np.arange(len(f))
frequency = f_range/TotalTime

n_kappa = int(np.floor(len(f_range)*1/max(frequency)))
fre_kappa = frequency[:n_kappa]
f_kappa = f[:n_kappa]

plt.xlabel('Frekvens [1/Y]')
plt.ylabel('Ampeltude [AU/Y]')
plt.grid()
plt.plot(fre_kappa, 2*abs(f_kappa)/N)
FrekI = []
for i in range (len(f_kappa)):
     if (2*abs(f_kappa[i])/N > 3e-6 ):
           print(f'angle = {np.angle(f[i])}, frekvens = {fre_kappa[i]}, ampeltude = {2*abs(f_kappa[i])/N }')
           FrekI.append(i)
plt.show()


plt.plot(TimeSample, NoiseSample ,'r')
signal =np.zeros(len(TimeSample))
for k in range(0,3):
      i = FrekI[k]
      signal += (2*abs(f_kappa[i])/N) *np.sin(2*np.pi*fre_kappa[i]*TimeSample + np.angle(f_kappa[i]) +np.pi/2)

plt.plot(TimeSample, signal,'gold')
plt.xlabel('Tid [Y]')
plt.ylabel('Fart [AU/Y]')
plt.show()