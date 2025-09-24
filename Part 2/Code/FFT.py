import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import  ast2000tools.constants as const

from scipy.ndimage import gaussian_filter


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
            q = np.exp(-1j *2 *np.pi/(N) *k)*odd[k]
            
            F[k] = (p + q)
            
            F[k + int(N/2)] = (p -q)
            
      return F

FileName = "SolarOrbitData.npz"
RawData = np.load(FileName)

TotalTime = RawData['time'][0]
N = 2**20

Data = RawData['v_noise'][:int(np.floor(430/TotalTime *len(RawData['v_raw'])))]
TotalTime = 430
v_pek = np.mean(Data)



NoiseData = Data -v_pek

TimeSample = np.linspace(0,TotalTime,N)
dt = TimeSample[1]- TimeSample[0]

index = np.linspace(0, len(NoiseData)-1,N, dtype=int)
NoiseSample = NoiseData[index] 


plt.plot(TimeSample, NoiseSample ,'r')
plt.plot(TimeSample, 0.00012*np.sin(2*np.pi*0.01627903*TimeSample) + 7.993e-6*np.sin(2*np.pi*0.3325*TimeSample) + 3.848e-6*np.sin(2*np.pi*0.05582*TimeSample),'gold')
plt.plot(TimeSample, 0.00012*np.sin(2*np.pi*0.01627903*TimeSample) + 7.993e-6*np.sin(2*np.pi*0.3325*TimeSample) ,'orange')
plt.plot(TimeSample, 0.00012*np.sin(2*np.pi*0.01627903*TimeSample) + 7.993e-6*np.sin(2*np.pi*0.3325*TimeSample) + 3.848e-6*np.sin(2*np.pi*0.05582*TimeSample)+2.341e-7*np.sin(2*np.pi*0.2223*TimeSample),'limegreen')

plt.show()


f = ditfft(NoiseSample,N)
f_range = np.arange(len(f))
frequency = f_range/TotalTime

n_kappa = int(np.floor(len(f_range)*1/max(frequency)))
fre_kappa = frequency[:n_kappa]
f_kappa = f[:n_kappa]

plt.plot(fre_kappa, 2*abs(f_kappa)/N)
plt.show()
