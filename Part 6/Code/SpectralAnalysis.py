# Regular imports
import numpy             as np
import matplotlib.pyplot as plt
from numba import njit
import pickle as pkl 

# AST imports
import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

"""
Our seed is [29578].
"""

k = const.k_B

@njit
def Sigma_l(l_0,k,t,m):
    return (l_0/(3*10**8)) * (((k*t)/m)**0.5)

@njit
def ChiSquare(fluxes, noises,l_0,m, T, F_min):

    O = np.zeros(3)
    X_min = 10e20

    lambda_steps = 1

    # Looping over all parameters
    for i_l in range(0,len(fluxes),lambda_steps):
        for i_t in range(len(T)):
            for i_f in range(len(F_min)):
                
                l_center = fluxes[i_l][0]
                t        = T[i_t]
                f_min    = F_min[i_f]

                X = 0

                # Looping over all fluxes
                for i in range(len(fluxes)):
                    l_i      = fluxes[i][0]
                    flux     = fluxes[i][1]
                    noise    = noises[i][1]

                    sigma_l = Sigma_l(l_0,k,t,m)
                    F_l = 1 + (f_min - 1)*np.exp(-((l_i-l_center)**2)/(2*sigma_l**2))
                    X += ((flux - F_l)/noise)**2

                if X < X_min:
                    X_min = X
                    O[0] = l_center
                    O[1] = t
                    O[2] = f_min

    return O

@njit
def ChiSquareSigma(fluxes, noises, S, F_min):

    O = np.zeros(3)
    X_min = 10e20

    # Looping over all parameters
    for i_l in range(len(fluxes)):
        for i_s in range(len(S)):
            for i_f in range(len(F_min)):
                
                l_center = fluxes[i_l][0]
                s        = S[i_s]
                f_min    = F_min[i_f]

                X = 0

                # Looping over all fluxes
                for i in range(len(fluxes)):
                    l_i      = fluxes[i][0]
                    flux     = fluxes[i][1]
                    noise    = noises[i][1]

                    F_l = 1 + (f_min - 1)*np.exp(-((l_i-l_center)**2)/(2*s**2))
                    X += ((flux - F_l)/noise)**2

                if X < X_min:
                    X_min = X
                    O[0] = l_center
                    O[1] = s
                    O[2] = f_min

    return O

njit
def ChiRange(fluxes, noises, L_0,M):

    T_min = 100; T_max = 500; N_T = 600//10
    F_min_min = 0.5; F_min_max = 1; N_F_min = int((F_min_max - F_min_min) * 500)

    T = np.linspace(T_min,T_max,N_T)
    F_min = np.linspace(F_min_min,F_min_max,N_F_min)

    O = np.zeros((len(L_0),3))

    for i in range(len(L_0)):

        i_l_min,i_l_max = FindLambda(L_0[i],len(fluxes)) # Note that these are indexes
    
        flux_range = fluxes[i_l_min:i_l_max]
        noise_range = noises[i_l_min:i_l_max]
        #S = np.linspace(Sigma_l(L_0[i],k,T_min,M[i]),Sigma_l(L_0[i],k,T_max,M[i]),N_T)

        print(f"Running optimization for lambda_0 : {L_0[i]:5}")
        O[i] = ChiSquare(flux_range,noise_range,L_0[i],M[i],T,F_min)

    return O

@njit
def FindLambda(l_0,length):
    d_l = (1/3) * l_0 * 10**(-4)

    l_max = l_0 + d_l
    l_min = l_0 - d_l

    i_max = int(((l_max - 600) / (2999.9998 - 600)) * length)
    i_min = int(((l_min - 600) / (2999.9998 - 600)) * length)

    return i_min,i_max

def ReadLists(flux_path, noise_path,step_length):
        
    num_steps = int(10000000/step_length)

    noises = np.zeros((num_steps,2))
    with open(noise_path,"r") as file:

        lines = file.readlines()
        i = 0
        for line in lines[::step_length]:

            line = line.strip(" ")
            words = line.split("     ")

            wavelength = float(words[0])
            noise = float(words[1])

            noises[i][0] = wavelength
            noises[i][1] = noise
            
            i += 1

    fluxes = np.zeros((num_steps,2))
    with open(flux_path,"r") as file:

        lines = file.readlines()
        i = 0
        for line in lines[::step_length]:

            line = line.strip(" ")
            words = line.split("     ")

            wavelength = float(words[0])
            flux = float(words[1])

            fluxes[i][0] = wavelength
            fluxes[i][1] = flux
            
            i += 1

    return fluxes,noises

if __name__ == "__main__":

    # Ast init
    seed = utils.get_seed('bmthune')
    mission = SpaceMission(seed)  
 
    # Files we will extract data from.
    noise_path = "sigma_noise.txt"
    flux_path = "spectrum_seed78_600nm_3000nm.txt"

    fluxes,noises = ReadLists(flux_path,noise_path,10)

    L_0 = np.array([ 
        632,
        690,
        760,
        720,
        820,
        940,
        1400,
        1600,
        1660,
        2200,
        2340,
        2870
        ])
    
    M   = np.array([ 
        32,
        32,
        32,
        18,
        18,
        18,
        44,
        44,
        16,
        16,
        28,
        44
        ]) * const.m_p
    
    print("\n---- Printing a few (lambda,F,sigma) triplets ----")
    N = 10
    Grain = 10001
    for i in range(Grain,N * (Grain+1),Grain):
        print(f"Lambda : {fluxes[i][0]:.5e}, F : {fluxes[i][1]:.5e}, sigma : {noises[i][1]}")
    print("---- ---------------------------------------- ----\n")

    # Plotting lambdas
    plt.plot(fluxes[:,0][:4000],fluxes[:,1][:4000],color = "limegreen")
    plt.axhline(1,color = "firebrick",linestyle = "--")
    plt.xlabel(r"Bølgelengde $\lambda$ [nm]")
    plt.ylabel(r"Relativ Fluks")
    plt.show()

    # Plotting sample spectra
    N = 1000
    X = np.linspace(0,1,N)
    x = 0.35
    def ColorMap(map,x):
        c = (x - np.min(x)) / (np.max(x) - np.min(x))
        return plt.get_cmap(map)(c)
    def GaussMap(c,mu,sigma,strength = 1):
        c_0 = c
        a = np.ones(len(c))
        for i in range(len(a)):
            x_i = i / len(a)
            y_i = ((2*np.pi*sigma**2)**0.5)*(1/((2*np.pi*sigma**2)**0.5)) * np.exp(-((x_i-mu)**2)/(2*sigma**2))
            a[i] = 1 - y_i
            if(a[i] < 0):
                a[i] = 0

        c_0[:,3] = a
        return c_0
    fig = plt.figure(figsize=(10, 8), dpi=80)
    ax = plt.axes()
    ax.bar(X,np.ones(N),width= 0.01,color = "k")
    colors = GaussMap(ColorMap("hsv",X),x,0.02)
    ax.bar(X,1,width= 0.01,color = colors) 
    ax.set_aspect(0.25)
    fig.tight_layout()
    plt.xticks([0,x,1],["",r"$\lambda_{center}$",""])
    plt.yticks([0,1],["",""])
    plt.show()

    print("Running Optimization")
    O = ChiRange(fluxes,noises,L_0,M)

    print("")

    for i in range(len(O)):
        print(f"Optimization for lambda_0 = {L_0[i]:5} : lambda = {O[i][0]:10.3f} nm, delta_lambda = {O[i][0] - L_0[i]:10.7f}, T = {O[i][1]:10.1f} K, F_min = {O[i][2]:20.3f}")
    print("")