import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import  ast2000tools.constants as const

#The Fast Fouier Transformation (FFT)
def ditfft(X, N):
      """The Fast Fouier Transformation, takes in a number of 
      steps and the number of datapoints sampled from the signal
      This is an algorithm to do Diskree Fourier Transformation"""
      if N == 1:
            return X


      

      even = ditfft(X[:N:2], int(N/2))
      odd = ditfft(X[1:N:2], int(N/2))

      F = np.zeros(N, dtype=np.complex128)
      for k in range (0, N//2):
            p = even[k]
            #This is the Cooley-Tukey algorithm, it splits the data points in even and odd and then splits the
            #sum (in a discre FT its a sum and not integralsign) Between even and odd, returning two new and equal
            #FT where one (q) is also multiplied with exp(-i 2 pi/n k), you can then express the original sum as a sum 
            #of theese two parts. For the integer over N/2 we have to a lot of algebra but end up that its even - odd (where odd is multiplied with the exp)
            q = np.exp(-1j * 2*np.pi/(N) * k) * odd[k] 
            F[k] = (p + q)
            F[k + int(N/2)] = (p -q)
      #Returns an array with the FT, aprox 0 in all places where frequency not found in data
      return F

#Loading Data containing the velocity of the star
FileName = "SolarOrbitData.npz"
RawData = np.load(FileName)
#Getting the total time for the duration 
TotalTime = RawData['time'][0]
#Choosing a number of steps, must be equal to a power of two for the algorithm to work
N = 2**20


Data = RawData['v_noise'][:int(np.floor(430/TotalTime *len(RawData['v_raw'])))] #Slicing the data so it contains an integer amount of periods
#If we dont slice like this the transform will give less accurate results

#This is the loongest timeperiod in this dataset that contains about an integer number of periods, we got this from reading of the plot
TotalTime = 430

#Fidning the pecular vellocity of the star and subtracting it from the array, the FT only works when centered around 0
v_pek = np.mean(Data)
Data = Data -v_pek

#Creating an array thats N loong and ecenly spread out integers between 0 and len Data
index = np.linspace(0, len(Data)-1,N, dtype=int)
#Sampling the Data with N evenly spread out datapoints
NoiseSample = Data[index] 


#Running the algorithm with inputs of the sampled data and N, returns the fouier transformed data f
f = ditfft(NoiseSample,N)

#first creating a range thats the len of f (N)
f_range = np.arange(len(f))
#The shortest frequency we can messure is 1/Totaltime, so this is the size of each frequency bin, then to find
#every value we can mesure we multiply it with the range above, this give freq = 1/Totaltime * range (1/totaltime is the fastest frequency it can messure)
frequency = f_range/TotalTime

#Halfing the size of both the frequency and f arrays because the second half is a mirrored coppy of the first half,
#This follows from the fact that the lowest frequency it can mesure is 1/2 totaltime not totaltime
n_kappa = int(np.floor(len(f_range)*1/max(frequency)))
fre_kappa = frequency[:n_kappa]
f_kappa = f[:n_kappa]

#Adding axies lables
plt.xlabel('Frekvens [1/Y]')
plt.ylabel('Ampeltude [AU/Y]')
#adding gridd
plt.grid()
plt.plot(fre_kappa, 2*abs(f_kappa)/N) #Plotting 2*|f|/N because this is the relation between the norm and the ampletude of the frequency
                                      #the plot therfore also shows the amplitude and the frequencies how this relation is derived I am not fully certain,
                                      #but is something i have read when recearching the FFT
FrekI = [] #creating an array to store the indexes of the outputs
for i in range (len(f_kappa)): #Here i loop over every ellement in the FT and check if there is an output over the noise
     if (2*abs(f_kappa[i])/N > 3.83e-6 ):
           #printing the amplitude frequency and Phase, exactly why the angle of the complex number is equal to the phase 
           #and the exact relation to the ampeltude are two things i have not fully grasped yet, but have read when reading about
           #the fft
           print(f'Phase = {np.angle(f[i])}, Frequency = {fre_kappa[i]}, Amplitude = {2*abs(f_kappa[i])/N }')
           FrekI.append(i)
plt.show()
#creating an array over the time with N steps
TimeSample = np.linspace(0,TotalTime,N)
#Ploting the signal, this can be ajusted over by changing the 'v_noise' to 'v_raw' to get the clean plot without noise
plt.plot(TimeSample, NoiseSample ,'r')
#creating an empty array of zeroes with len N (timesteps)
signal =np.zeros(len(TimeSample))
for i in FrekI:
      #Creating the singal that contains all the values that come from the FT, both phase amplitude and frequency
      signal += (2*abs(f_kappa[i])/N) *np.sin(2*np.pi*fre_kappa[i]*TimeSample + np.angle(f_kappa[i]) +np.pi/2)

plt.plot(TimeSample, signal,'gold')
#Naming the axies
plt.xlabel('Tid [Y]')
plt.ylabel('Fart [AU/Y]')
#Showing
plt.show()