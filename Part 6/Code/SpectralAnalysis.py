# BRUKER IKKE KODEMAL
# Skrevet av Bastian Eggum Huuse og Bendik Thune

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
Our seed is [29578], meaning we need the file
spectrum_seed78_600nm_3000nm.txt, along with
sigma_noise.txt
"""

# Saving the boltzmann-constant
k = const.k_B

@njit
def Sigma_l(l_0,k,t,m):
    """
    Function that calculates the standard deviation for a line profile around the 
    wavelength l_0

    parameters:
    l_0 : float | The wavelength we build our line profile around
    k   : float | The boltzmann constant
    t   : float | The assumed temperature of the gas
    m   : float | The mass of the particle that absorbs l_0

    returns:
    float | standard deviation for the line profile
    """

    return (l_0/(3*10**8)) * (((k*t)/m)**0.5)

@njit
def ChiSquare(fluxes, noises,l_0,m, T, F_min):
    """
    
    Function that performs a Chi Square optimization, fitting a flux curve to the 
    measured flux points. Note that fluxes and noises are sliced before being passed into this function,
    so they don't contain all flux and noise values.

    parameters:
    fluxes : Array(float) | Array containing (wavelength,measured flux value) pairs
    noises : Array(float) | Array containing (wavelength,measured noise value) pairs
    l_0    : float        | The wavelength we wish to construct a line profile for (before doppler-shift)
    m      : float        | The mass of the particle that absorbs l_0
    T      : Array(float) | Array containing all possible Temperature values for the gas
    F_min  : Array(float) | Array containing all possible F_min values for the line profile
    
    returns:
    Array(float) | A triplet containing the parameters of the optimized line profile function (l_c,t,f_min)
    """

    # Defining our O array, containing our optimized parameters, and a variable
    # containing the lowest Chi-value. We initialize this as a huge number
    O = np.zeros(3)
    X_min = 10e20

    # Looping over all parameters
    for i_l in range(0,len(fluxes)):
        for i_t in range(len(T)):
            for i_f in range(len(F_min)):
                
                # Grabbing the current set of parameters
                l_center = fluxes[i_l][0]
                t        = T[i_t]
                f_min    = F_min[i_f]

                # Initializing a Chi value to sum over
                X = 0
                sigma_l = Sigma_l(l_0,k,t,m)

                # Looping over all fluxes
                for i in range(len(fluxes)):

                    # Grabbing current Wavelength, flux, and noise
                    l_i      = fluxes[i][0]
                    flux     = fluxes[i][1]
                    noise    = noises[i][1]

<<<<<<< HEAD
=======
                    # Calculating sigma
                    sigma_l = Sigma_l(l_0,k,t,m)
                    # Calculating line profile flux
>>>>>>> e34fba0d0bc5e3ef62f6ba2d0c43adafae0ea113
                    F_l = 1 + (f_min - 1)*np.exp(-((l_i-l_center)**2)/(2*sigma_l**2))
                    # Adding to Chi
                    X += ((flux - F_l)/noise)**2

                # If the current Chi is the smalles one we've found, replace the parameters
                if X < X_min:
                    X_min = X
                    O[0] = l_center
                    O[1] = t
                    O[2] = f_min

    return O

njit
def ChiRange(fluxes, noises, L_0,M):

    """
    Function that runs a chi-square optimization for many l_0 values.

    parameters:
    fluxes : Array(float) | Array containing all (wavelength,measured flux value) pairs
    noises : Array(float) | Array containing all (wavelength,measured noise value) pairs
    L_0    : Array(float) | Array containing all l_0 values we want to find line profiles for
    M      : Array(float) | Array containing the masses of the particles that absorb the l_0 values

    returns:
    Array(float) | Array containing all sets of optimized parameters. (l_c,t,f_min) 
    """

    # Defining boundaries for T and F_min values.
    T_min     = 100; T_max     = 500; N_T     = 600//10
    F_min_min = 0.5; F_min_max = 1.0; N_F_min = int((F_min_max - F_min_min) * 500)

    # Creating arrays of T and F_min values
    T     = np.linspace(T_min,T_max,N_T)
    F_min = np.linspace(F_min_min,F_min_max,N_F_min)

    # Creating array to contain optimized parameters
    O = np.zeros((len(L_0),3))

    # Looping over all lambda_0 values
    for i in range(len(L_0)):

        # Computing max and min lambda index (within fluxes)
        i_l_min,i_l_max = FindLambda(L_0[i],len(fluxes)) # Note that these are indexes
    
        # Slicing fluxes and noises
        flux_range = fluxes[i_l_min:i_l_max]
        noise_range = noises[i_l_min:i_l_max]

        # Running optimization for current l_0 value
        print(f"Running optimization for lambda_0 : {L_0[i]:5}")
        O[i] = ChiSquare(flux_range,noise_range,L_0[i],M[i],T,F_min)

    return O

@njit
def FindLambda(l_0,length):

    """
    Short function that computes the indexes of the maximum and minimum lambda values

    parameters:
    l_0    : float | The lambda value we want to create a boundary around.
    length : float | The total lenght of the fluxes array

    returns:
    int | index corresponding to the minimum wavelength we want to search at
    int | index corresponding to the maximum wavelength we want to search at
    """

    d_l = (1/3) * l_0 * 10**(-4)

    l_max = l_0 + d_l
    l_min = l_0 - d_l

    i_max = int(((l_max - 600) / (2999.9998 - 600)) * length)
    i_min = int(((l_min - 600) / (2999.9998 - 600)) * length)

    return i_min,i_max

def PlotLineprofile(fluxes,l_0,m,O):

    """
    Function that plots the graphs of the measured flux values, along with the ones
    from the line profiles
    """

    i_min,i_max = FindLambda(l_0,len(fluxes))
    trueFluxes = fluxes[i_min:i_max]

    sigma_l = Sigma_l(l_0,k,O[1],m)
    F_l = 1 + (O[2] - 1)*np.exp(-((trueFluxes[:,0]-O[0])**2)/(2*sigma_l**2))

    ax = plt.axes()
    ax.plot(trueFluxes[:,0],trueFluxes[:,1])
    ax.plot(trueFluxes[:,0],F_l)

    plt.title(f"Linjeprofil for lambda_0 = {l_0} nm")
    plt.xlabel(r"Bølgelengde $\lambda$")
    plt.ylabel("Relativ fluks")
    plt.show()

if __name__ == "__main__":

    # Ast init
    seed = utils.get_seed('bmthune')
    mission = SpaceMission(seed)  
 
    # Files we will extract data from.
    noise_path = "sigma_noise.txt"
    flux_path = "spectrum_seed78_600nm_3000nm.txt"

    """ Reading files """

    fluxes = np.loadtxt(flux_path)
    noises = np.loadtxt(noise_path)
    step_length = 1
    fluxes = fluxes[::step_length]
    noises = noises[::step_length]

    """ Defining arrays for lambda_0 values and corresponding masses """

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
    

    """ Printing and plotting some preliminary info"""

    print("\n---- Printing a few (lambda,F,sigma) triplets ----")
    N = 10
    Grain = 100010
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

    """ Actually running the optimization """
    print("Running Optimization")
    O = ChiRange(fluxes,noises,L_0,M)
    print("")


    """ Printing and plotting results """
    for i in range(len(O)):
        print(f"Optimization for lambda_0 = {L_0[i]:5} : lambda = {O[i][0]:10.3f} nm, delta_lambda = {O[i][0] - L_0[i]:10.7f}, T = {O[i][1]:10.1f} K, F_min = {O[i][2]:20.3f}")
    print("")
    for i in range(len(L_0)):
        PlotLineprofile(fluxes,L_0[i],M[i],O[i])

"""
Output:

---- Printing a few (lambda,F,sigma) triplets ----
Lambda : 6.24002e+02, F : 9.84610e-01, sigma : 0.050055427
Lambda : 6.48005e+02, F : 1.01006e+00, sigma : 0.050179341
Lambda : 6.72007e+02, F : 9.61666e-01, sigma : 0.050268934
Lambda : 6.96010e+02, F : 9.65090e-01, sigma : 0.050221439
Lambda : 7.20012e+02, F : 9.57286e-01, sigma : 0.050000413
Lambda : 7.44014e+02, F : 9.86879e-01, sigma : 0.050332957
Lambda : 7.68017e+02, F : 9.86082e-01, sigma : 0.050627867
Lambda : 7.92019e+02, F : 9.45029e-01, sigma : 0.05071695
Lambda : 8.16022e+02, F : 1.02570e+00, sigma : 0.050497633
Lambda : 8.40024e+02, F : 9.08074e-01, sigma : 0.050001649
---- ---------------------------------------- ----

Running Optimization
Running optimization for lambda_0 :   632
Running optimization for lambda_0 :   690
Running optimization for lambda_0 :   760
Running optimization for lambda_0 :   720
Running optimization for lambda_0 :   820
Running optimization for lambda_0 :   940
Running optimization for lambda_0 :  1400
Running optimization for lambda_0 :  1600
Running optimization for lambda_0 :  1660
Running optimization for lambda_0 :  2200
Running optimization for lambda_0 :  2340
Running optimization for lambda_0 :  2870

Optimization for lambda_0 =   632 : lambda =    632.005 nm, delta_lambda =  0.0047200, T =      106.8 K, F_min =                0.851
Optimization for lambda_0 =   690 : lambda =    690.007 nm, delta_lambda =  0.0074400, T =      500.0 K, F_min =                0.765
Optimization for lambda_0 =   760 : lambda =    760.022 nm, delta_lambda =  0.0216800, T =      100.0 K, F_min =                0.863
Optimization for lambda_0 =   720 : lambda =    720.020 nm, delta_lambda =  0.0199200, T =      100.0 K, F_min =                0.940
Optimization for lambda_0 =   820 : lambda =    820.003 nm, delta_lambda =  0.0032000, T =      100.0 K, F_min =                0.789
Optimization for lambda_0 =   940 : lambda =    939.987 nm, delta_lambda = -0.0128800, T =      100.0 K, F_min =                0.894
Optimization for lambda_0 =  1400 : lambda =   1399.990 nm, delta_lambda = -0.0102000, T =      100.0 K, F_min =                0.871
Optimization for lambda_0 =  1600 : lambda =   1599.958 nm, delta_lambda = -0.0422000, T =      276.3 K, F_min =                0.805
Optimization for lambda_0 =  1660 : lambda =   1660.029 nm, delta_lambda =  0.0286000, T =      500.0 K, F_min =                0.946
Optimization for lambda_0 =  2200 : lambda =   2199.992 nm, delta_lambda = -0.0078000, T =      161.0 K, F_min =                0.894
Optimization for lambda_0 =  2340 : lambda =   2340.040 nm, delta_lambda =  0.0396000, T =      100.0 K, F_min =                0.942
Optimization for lambda_0 =  2870 : lambda =   2869.924 nm, delta_lambda = -0.0757000, T =      181.4 K, F_min =                0.801
"""