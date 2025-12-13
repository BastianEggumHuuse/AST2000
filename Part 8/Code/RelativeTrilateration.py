# Ikke brukt kodemal!!
# Skrevet av Bastian Eggum Huuse og Bendik Thune

import sys
import numpy as np
from numba import njit
import matplotlib.pyplot as plt

import ast2000tools.constants as const
import ast2000tools.utils     as utils

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

def TrilaterationAlgorithm(r_p,Positions, Distances):

    # Defining star-distance and the range of angles
    Grain = 15
    Range = np.linspace(0,2*np.pi,2**(Grain))

    # Running the algorithm to find the angle from the sun at which our position is at.
    Angle = BinaryLeastSquares(r_p,Range,Positions,Distances)

    # We want to return the position, instead of just the angle
    return (r_p * np.cos(Angle),r_p * np.sin(Angle))