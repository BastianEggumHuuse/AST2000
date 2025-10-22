# Ikke brukt kodemal!!
# Skrevet av Bastian Eggum Huuse og Bendik Thune

import sys
import numpy as np
from numba import njit
import matplotlib.pyplot as plt

import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

from GeneralizedLaunch import NumericalOrbitFunction,main

@njit
def TupleDiff(Tuple1,Tuple2):

    """
    Support function that returns the difference between two tuples containing float values
    (Note that both tuple functions also work with arrays)

    Parameters:
    Tuple1 : tuple | The first tuple
    Tuple2 : tuple | The second duple

    Returns:
    tuple | The difference between the two tuples
    """

    Coord0 = Tuple1[0] - Tuple2[0]
    Coord1 = Tuple1[1] - Tuple2[1]
    return Coord0,Coord1

@njit
def TupleNorm(Tuple):

    """
    Support function that returns the norm of a tuple containing two float values

    Parameters:
    Tuple : tuple | The tuple

    Returns:
    float | The norm of the tuple
    """

    return (Tuple[0]**2 + Tuple[1]**2)**(1/2)

@njit
def PlanetDifference(Pos,PlanetPositions,TruePlanetDistances):

    """
    Method that computes the difference between the distances from the planets to a given point, and what those distances should be.
    
    Parameters:
    Pos                 : Array(float) | The current position we want to check
    PlanetPositions     : Array(float) | List of planet positions
    TruePlanetDistances : Array(float) | List of expected distances

    returns: float | The difference between computed distances and expected distances (the sum of these)
    """

    PlanetDeviations = np.zeros((len(PlanetPositions),2))

    # Finding the deviation vector between the given point and all planets
    for i in range(len(PlanetDeviations)):
        Deviation = TupleDiff(Pos,PlanetPositions[i])
        PlanetDeviations[i][0] = Deviation[0]
        PlanetDeviations[i][1] = Deviation[1]

    TotalDifference = 0

    # Computing the sum of the difference between each computed distance, and each expected distance
    for i in range(len(PlanetDeviations)):
        TotalDifference += abs(TupleNorm(PlanetDeviations[i]) - TruePlanetDistances[i])

    return TotalDifference

@njit
def BinaryLeastSquares(StarDistance,Range,PlanetPositions,TruePlanetDistances):

    """
    Algorithm that finds the point where expected distances match computed distances (closest). This would then be the point in space where the rocket is.
    This point is found through a combination binary search and least squares approach.

    We create a circle with center in the star, and radius equal to the distance the rocket is from the sun.
    We place a point on each side of this circle, and find the point closest to the expected position through least squares
    We then select the half of the circle the closest point is on, and then call the method recursively on that part of the circle.

    Parameters:
    StarDistance        : float        | The distance between the rocket and the star
    Range               : float        | The range of angles we are currently checking (Begins as 0-2pi)
    PlanetPositions     : Array(float) | The positions of all the planets
    TruePlanetDistances : Array(float) | The expected distances of each planet 
    
    returns: float | The angle from the star the position of the rocket is located (in radians)
    """

    # If the range of angles has a length of 1, we simply return that angle (We've found the closest angle)
    if len(Range) == 1:
        return Range[0]

    # Splitting Range in two
    HalfLength = int(len(Range)/2)
    Range1 = np.zeros(HalfLength)
    Range2 = np.zeros(HalfLength)

    # Looping through and duplicating the half ranges
    for i in np.arange(HalfLength*2):
        
        if i < HalfLength:
            Range1[i] = Range[i]
        else:
            Range2[i - HalfLength] = Range[i]

    # Finding the Positions we want to check
    HalfIndex = int(HalfLength/2)
    Position1 = (StarDistance * np.cos(Range1[HalfIndex]),StarDistance * np.sin(Range1[HalfIndex]))
    Position2 = (StarDistance * np.cos(Range2[HalfIndex]),StarDistance * np.sin(Range2[HalfIndex]))

    # Finding the least squares difference for both of these points
    Diff1 = PlanetDifference(Position1,PlanetPositions,TruePlanetDistances)
    Diff2 = PlanetDifference(Position2,PlanetPositions,TruePlanetDistances)

    # Recursing with the smallest one
    if Diff1 <= Diff2:
        Angle = BinaryLeastSquares(StarDistance,Range1,PlanetPositions,TruePlanetDistances)
    else:
        Angle = BinaryLeastSquares(StarDistance,Range2,PlanetPositions,TruePlanetDistances)

    return Angle

def TrilaterationAlgorithm(t,Distances):
    
    """
    Main function that actually calls the algorithm.

    Parameters:
    t         : float        | Time of launch (in years)
    Distances : Array(float) | Array containing the Expected distances of all planets (and the star)

    returns: Array(float) | The position of the rocket after launch (in AU)
    """

    # Instantiating our Position function
    FilePath = "NumericalOrbitData.npz"
    PlanetPositionFunction = NumericalOrbitFunction(FilePath)

    # Defining star-distance and the range of angles
    StarDistance = Distances[-1]
    Grain = 10
    Range = np.linspace(0,2*np.pi,2**(Grain))
    
    # Defining an array of distances that doesn't include the sun
    PlanetDistances = np.zeros((len(Distances)-1))
    for i in range(len(PlanetDistances)):
        PlanetDistances[i] = Distances[i]

    # Defining an array of positions, one for each planet.
    PlanetPositions = np.zeros((len(PlanetDistances),2))
    for p in range(len(PlanetPositions)):
        PlanetPositions[i] = PlanetPositionFunction(t,p)

    # Running the algorithm to find the angle from the sun at which our position is at.
    Angle = BinaryLeastSquares(StarDistance,Range,PlanetPositions,PlanetDistances)

    # We want to return the position, instead of just the angle
    return (StarDistance * np.cos(Angle),StarDistance * np.sin(Angle))

if __name__ ==  "__main__":

    # Ast init
    seed = utils.get_seed('bmthune')
    mission = SpaceMission(seed)

    t_0 = 1

    # Running Launch
    t_1 = main(mission,t_0)

    # Collecting Distances
    Distances = mission.measure_distances()

    # Running Algorithm
    Position = TrilaterationAlgorithm(t_1,Distances)

    print(f"\nPosition after launch (Computed from distances): {Position} AU")




    # Plotting
    # Initializing plotting
    fig, ax = plt.subplots()

    FilePath = "NumericalOrbitData.npz"
    PlanetPositionFunction = NumericalOrbitFunction(FilePath)

    # Defining an array of distances that doesn't include the sun
    PlanetDistances = np.zeros((len(Distances)-1))
    for i in range(len(PlanetDistances)):
        PlanetDistances[i] = Distances[i]

    # Defining an array of positions, one for each planet.
    PlanetPositions = np.zeros((len(PlanetDistances),2))
    for p in range(len(PlanetPositions)):
        PlanetPositions[p] = PlanetPositionFunction(t_1,p)

    # Adding the star circle
    star = plt.Circle((0, 0), Distances[-1], color = 'gold',fill = False)
    ax.add_patch(star)

    for i in range(len(PlanetDistances)):
        # Adding the star (not to scale)
        planet = plt.Circle(PlanetPositions[i], PlanetDistances[i],linestyle = "-",fill = False)
        ax.add_patch(planet)

    ax.plot(0,0,"*",color = "gold")

    # Making axes equal, and showing plot
    plt.xlabel("Posisjon på x-aksen [AU]")
    plt.ylabel("Posisjon på y-aksen [AU]")
    plt.title("Skjæring mellom sirkler med radius lik distanse til rakett")
    plt.axis('equal')
    plt.show()

"""
Code that runs this program: python Trilateration.py

This program also outputs a plot, which has been turned in as Trilateration_plot.py
Output:

Initializing motor...
Calculated Force per motor : 1.26686e-10, Calculated Fuel Consumption per motor : 2.95044e-14
Calculated Force           : 2.11143e+06, Calculated Fuel Consumption      : 4.91740e+02

Generalized Position in solar system frame at t = 1: [x : -1.302 AU, y : 2.42e+00 AU]
Generalized Velocity in solar system frame at t = 1: [x : -6.414 AU/Y, y : -0.728 AU/Y]

Specialized Position in solar system frame at t = 1: [x : 2.814 AU, y : 6.98e-05 AU]
Specialized Velocity in solar system frame at t = 1: [x : 2.480 AU/Y, y : 5.856 AU/Y]

Rocket was moved down by 6524.35 m to stand on planet surface.
New launch parameters set.
Launch completed, reached escape velocity in 379.294 s.
Your spacecraft position was satisfyingly calculated. Well done!
*** Achievement unlocked: No free launch! ***

Position after launch (Computed from distances): (-1.3002171920399341, 2.420787581930546) AU
"""