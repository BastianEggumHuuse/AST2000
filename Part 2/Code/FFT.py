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
            
      
      e = ditfft(X[:N:2], int(N/2))
      o = ditfft(X[1:N:2], int(N/2))

      for k in range (0, N//2):
            p = e[k]
            q = np.exp(-1j *2 *np.pi/(N) *k)*o[k]
            
            X[k] = (p + q)
            
            X[k + int(N/2)] = (p -q)
            
      return X

FileName = "SolarOrbitData.npz"
Data = np.load(FileName)

NoiseData = Data["v_noise"]
RawData   = Data["v_raw"]
TotalTime = Data["time"][0]

v_pec = np.mean(RawData)
RawData -= v_pec
N = 2**15

Plateau = np.ones(len(NoiseData))
ratio = 0.05
Indexes = np.arange(int(np.floor(ratio * len(Plateau))))
LerpersIn = np.sin(np.linspace(0,np.pi/2,len(Indexes)))
Plateau[Indexes] *= LerpersIn
Plateau[len(Plateau) - Indexes - 1] = LerpersIn

plt.plot(np.linspace(0,1,len(NoiseData)),Plateau)

TimeSample = np.linspace(0,TotalTime,len(NoiseData))[:N]
NoiseSample = RawData[:N].astype(np.complex128)
f = ditfft(NoiseSample,N)

f_range = np.arange(len(f))
frequency = f_range/TimeSample[-1]

#plt.plot(frequency,abs(f))
plt.show()