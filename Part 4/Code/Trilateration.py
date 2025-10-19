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

    Coord0 = Tuple1[0] - Tuple2[0]
    Coord1 = Tuple1[1] - Tuple2[1]
    return Coord0,Coord1

@njit
def TupleNorm(Tuple):

    return (Tuple[0]**2 + Tuple[1]**2)**(1/2)

@njit
def PlanetDifference(Pos,PlanetPositions,TruePlanetDistances):

    PlanetDeviations = np.zeros((len(PlanetPositions),2))

    for i in range(len(PlanetDeviations)):
        Deviation = TupleDiff(Pos,PlanetPositions[i])
        PlanetDeviations[i][0] = Deviation[0]
        PlanetDeviations[i][1] = Deviation[1]

    TotalDifference = 0

    for i in range(len(PlanetDeviations)):
        TotalDifference += abs(TupleNorm(PlanetDeviations[i]) - TruePlanetDistances[i])

    return TotalDifference

@njit
def BinaryLeastSquares(StarDistance,Range,PlanetPositions,TruePlanetDistances):

    if len(Range) == 1:
        return Range[0]

    # Splitting Range in two
    HalfLength = int(len(Range)/2)
    Range1 = np.zeros(HalfLength)
    Range2 = np.zeros(HalfLength)

    for i in np.arange(HalfLength*2):
        
        if i < HalfLength:
            Range1[i] = Range[i]
        else:
            Range2[i - HalfLength] = Range[i]

    # Finding the distances
    HalfIndex = int(HalfLength/2)
    Position1 = (StarDistance * np.cos(Range1[HalfIndex]),StarDistance * np.sin(Range1[HalfIndex]))
    Position2 = (StarDistance * np.cos(Range2[HalfIndex]),StarDistance * np.sin(Range2[HalfIndex]))

    Diff1 = PlanetDifference(Position1,PlanetPositions,TruePlanetDistances)
    Diff2 = PlanetDifference(Position2,PlanetPositions,TruePlanetDistances)

    if Diff1 <= Diff2:
        Angle = BinaryLeastSquares(StarDistance,Range1,PlanetPositions,TruePlanetDistances)
    else:
        Angle = BinaryLeastSquares(StarDistance,Range2,PlanetPositions,TruePlanetDistances)

    return Angle

def TrilaterationAlgorithm(t,Distances):
    
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

    # Run Launch
    t_1 = main(mission,t_0)

    Distances = mission.measure_distances()

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
