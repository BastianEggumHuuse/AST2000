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

def PlanetDifference(Pos,FirstPlanetPos,TrueFirstPlanetDist,SecondPlanetPos,TrueSecondPlanetDist):

    FirstPlanetDeviation = TupleDiff(FirstPlanetPos,Pos)
    SecondPlanetDeviation = TupleDiff(FirstPlanetPos,Pos)

    FirstPlanetDifference = abs(TupleNorm(FirstPlanetDeviation) - TrueFirstPlanetDist)
    SecondPlanetDifference = abs(TupleNorm(SecondPlanetDeviation) - TrueSecondPlanetDist)

    TotalDifference = FirstPlanetDifference + SecondPlanetDifference

    return TotalDifference

def BinaryLeastSquares(StarDistance,Range,FirstPlanetPos,SecondPlanetPos,TrueFirstPlanetDist,TrueSecondPlanetDist):

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

    Diff1 = PlanetDifference(Position1,FirstPlanetPos,TrueFirstPlanetDist,SecondPlanetPos,TrueSecondPlanetDist)
    Diff2 = PlanetDifference(Position2,FirstPlanetPos,TrueFirstPlanetDist,SecondPlanetPos,TrueSecondPlanetDist)

    if Diff1 <= Diff2:
        Angle = BinaryLeastSquares(StarDistance,Range1,FirstPlanetPos,SecondPlanetPos,TrueFirstPlanetDist,TrueSecondPlanetDist)
    else:
        Angle = BinaryLeastSquares(StarDistance,Range2,FirstPlanetPos,SecondPlanetPos,TrueFirstPlanetDist,TrueSecondPlanetDist)

    return Angle

def Main(t,Distances):
    
    FilePath = "NumericalOrbitData.npz"
    PlanetPositionFunction = NumericalOrbitFunction(FilePath)

    StarDistance = Distances[-1]

    FirstPlanetPos = PlanetPositionFunction(t,0)
    SecondPlanetPos = PlanetPositionFunction(t,1)

    TrueFirstPlanetDist = Distances[0]
    TrueSecondPlanetDist = Distances[1]

    Range = np.linspace(0,2*np.pi,2**(10))

    Angle = BinaryLeastSquares(StarDistance,Range,FirstPlanetPos,SecondPlanetPos,TrueFirstPlanetDist,TrueSecondPlanetDist)

    return (StarDistance * np.cos(Angle),StarDistance * np.sin(Angle))