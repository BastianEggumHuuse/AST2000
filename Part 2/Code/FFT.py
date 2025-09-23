import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import  ast2000tools.constants as const

seed = utils.get_seed('bmthune')
system = SolarSystem(seed)

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

FileName = "SolarOrbitData.npz"
Data = np.load(FileName)

NoiseData = Data["v_noise"]
RawData   = Data["v_raw"]
TotalTime = Data["time"][0]
N = 2
while N*2 < len(NoiseData):
      N *= 2
      
TimeSample = np.linspace(0,TotalTime,len(NoiseData))[:N]
NoiseSample = RawData[:N]
f = ditfft(NoiseSample,N)

f_range = np.arange(len(f))
frequency = f_range/TimeSample[-1]

plt.plot(frequency,abs(f))
plt.show()