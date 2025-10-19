# Ikke brukt kodemal!!
# Skrevet av Bastian Eggum Huuse

import sys
import numpy as np
from numba import njit

import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

from GeneralizedLaunch import NumericalOrbitFunction

def TupleDiff(Tuple1,Tuple2):

    Coord0 = Tuple1[0] - Tuple2[0]
    Coord1 = Tuple1[1] - Tuple2[1]
    return (Coord0,Coord1)

def TupleNorm(Tuple):

    return (Tuple[0]**2 + Tuple[1]**2)**(1/2)

def PlanetDifference(Pos,PlanetPositions,TruePlanetDistances):

    PlanetDeviations = np.zeros((len(PlanetPositions,2)))

    for i in range(len(PlanetDeviations)):
        PlanetDeviations[i] = TupleDiff(PlanetPositions,Pos)

    TotalDifference = 0

    for i in range(len(PlanetDeviations)):
        TotalDifference += abs(TupleNorm(PlanetDeviations[i]) - TruePlanetDistances[i])

    return TotalDifference

def BinaryLeastSquares(StarDistance,Range,PlanetPositions,TruePlanetDistances):

    if len(Range) == 1:
        return Range[0]

    # Splitting Range in two
    HalfLength = len(Range)/2
    Range1 = np.zeros(HalfLength)
    Range2 = np.zeros(HalfLength)

    for i in np.arange(HalfLength*2):
        
        if i < HalfLength:
            Range1[i] = Range[i]
        else:
            Range2[i - HalfLength] = Range[i]

    # Finding the distances
    Position1 = (StarDistance * np.cos(Range1[HalfLength/2]),StarDistance * np.sin(Range1[HalfLength/2]))
    Position2 = (StarDistance * np.cos(Range2[HalfLength/2]),StarDistance * np.sin(Range2[HalfLength/2]))

    Diff1 = PlanetDifference(Position1,PlanetPositions,TruePlanetDistances)
    Diff2 = PlanetDifference(Position2,PlanetPositions,TruePlanetDistances)

    if Diff1 <= Diff2:
        Angle = BinaryLeastSquares(StarDistance,Range1,PlanetPositions,TruePlanetDistances)
    else:
        Angle = BinaryLeastSquares(StarDistance,Range2,PlanetPositions,TruePlanetDistances)

    return Angle

def Main(t,Distances):
    
    FilePath = "NumericalOrbitData.npz"
    PlanetPositionFunction = NumericalOrbitFunction(FilePath)

    # Defining star-distance and the range of angles
    StarDistance = Distances[-1]
    Range = np.linspace(0,2*np.pi,2**(10))
    
    # Defining an array of distances that doesn't include the sun
    PlanetDistances = np.zeros((len(Distances)-1))
    for i in range(len(Distances)):
        PlanetDistances[i] = Distances[i]

    # Defining an array of positions, one for each planet.
    PlanetPositions = np.zeros((len(PlanetDistances),2))
    for p in range(len(PlanetPositions)):
        PlanetPositions[i] = PlanetPositionFunction(t,p)

    # Running the algorithm to find the angle from the sun at which our position is at.
    Angle = BinaryLeastSquares(StarDistance,Range,Distances)

    # We want to return the position, instead of just the angle
    return (StarDistance * np.cos(Angle),StarDistance * np.sin(Angle))

if __name__ ==  "__main__":

    # Ast init
    seed = utils.get_seed('bmthune')
    mission = SpaceMission(seed)

    # Run Launch
    
    Distances = mission.measure_distances()

    