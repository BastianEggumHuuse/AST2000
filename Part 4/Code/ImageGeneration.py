# Ikke brukt kodemal!!
# Skrevet av Bastian Eggum Huuse

import sys
import numpy as np
from numba import njit
from PIL import Image

import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

@njit
def GetImageDimentions(Pixels):
    
    Height = len(Pixels)
    Width = len(Pixels[0])

    return Height,Width

@njit
def MaxFunction(Alpha):
    
    Value = (2*np.sin(Alpha/2))/(1+np.cos(Alpha/2))
    return Value

@njit
def MaxMin(FOV):
    
    Value = MaxFunction(FOV)

    Max = Value
    Min = -Value

    return Min,Max

@njit
def Rho(x,y):
    return (x**2 + y**2)**(1/2)

@njit
def Beta(x,y):
    rho = Rho(x,y)
    return 2* np.arctan(rho/2)

@njit
def Phi(phi0,theta0,x,y):

    rho = Rho(x,y)
    beta = Beta(x,y)

    return phi0 + np.arctan((x * np.sin(beta))/(rho*np.sin(theta0) * np.cos(beta) - y * np.cos(theta0)*np.sin(beta)))

@njit
def Theta(phi0,theta0,x,y):

    rho = Rho(x,y)
    beta = Beta(x,y)

    return theta0 - np.arcsin(np.cos(beta)*np.cos(theta0) + (y/rho)*np.sin(beta)*np.sin(theta0))

@njit
def GenerateCoordinateGrid(RangeX,RangeY,PixelWidth,PixelHeight):

    ValuesX = np.linspace(RangeX[0],RangeX[1],PixelWidth)
    ValuesY = np.linspace(RangeY[0],RangeY[1],PixelHeight)
    # The image gives us height as the first dimention, but we want width as the first dimention
    CoordinateGrid = np.zeros((PixelWidth,PixelHeight,2))

    for ix in range(len(ValuesX)):
        for iy in range(len(ValuesY)):
            CoordinateGrid[ix][iy] = np.array([ValuesX[ix],ValuesY[iy]])

    return CoordinateGrid

@njit
def GenerateAngleGrid(CoordinateGrid,phi0,theta0,Wrap = False):
    
    AngleGrid = np.zeros((len(CoordinateGrid),len(CoordinateGrid[0]),2))

    for ix in range(len(CoordinateGrid)):
        for iy in range(len(CoordinateGrid[ix])):
            
            x = CoordinateGrid[ix][iy][0]
            y = CoordinateGrid[ix][iy][1]

            theta = Theta(phi0,theta0,x,y)
            phi = Phi(phi0,theta0,x,y)

            if(Wrap):
                if theta > np.pi:
                    theta = np.pi
                while theta < 0:
                    theta = 0

            AngleGrid[ix][iy] = np.array([theta,phi])

    return AngleGrid

def GenerateImage(SkyData,AngleGrid):
    
    ImageGrid = np.zeros((AngleGrid.shape[1],AngleGrid.shape[0],3))

    for ix in range(len(AngleGrid)):
        for iy in range(len(AngleGrid[0])):
            PixelIndex = SpaceMission.get_sky_image_pixel(AngleGrid[ix][iy][0],AngleGrid[ix][iy][1])
            ImageGrid[iy][ix] = SkyData[PixelIndex][2:]

    return ImageGrid.astype('uint8')

def GenerateImageRange(Skydata,CoordinateGrid,theta0):

    phi0 = np.arange(0,360,1).astype("float32")

    Images = np.zeros((len(phi0),CoordinateGrid.shape[1],CoordinateGrid.shape[0],3),dtype="uint8")

    for i in range(len(phi0)):

        phi0[i] = np.deg2rad(phi0[i])

        AngleGrid = GenerateAngleGrid(CoordinateGrid,phi0[i],theta0)
        ImageGrid = GenerateImage(Skydata,AngleGrid)
        Images[i] = ImageGrid

    return Images

def Main():

    SkyData = np.load("himmelkule.npy")

    InputImage = Image.open("sample0000.png")
    FOV_PHI = np.deg2rad(70)
    FOV_THETA = np.deg2rad(70)

    # Processing image
    InputImage = InputImage
    Pixels = np.array(InputImage)
    PixelHeight, PixelWidth = GetImageDimentions(Pixels)

    # Generating Ranges
    RangeX = np.array(MaxMin(FOV_PHI))
    RangeY = -np.array(MaxMin(FOV_THETA))

    # Generating grid of x and y values:
    CoordinateGrid = GenerateCoordinateGrid(RangeX,RangeY,PixelWidth,PixelHeight)
    AngleGrid = GenerateAngleGrid(CoordinateGrid,phi0 = 0,theta0 = np.pi/2)
    # Base parameters : phi0 = 0 | theta0 = np.pi/2
    ImageGrid = GenerateImage(SkyData,AngleGrid)

    # Reproducing sample2000.png
    OutputImage = Image.fromarray(ImageGrid)
    OutputImage.save("Output0000.png") # Make new png

    # running for all phi0 between 0 and 2pi

    ImageGrids = GenerateImageRange(SkyData,CoordinateGrid,theta0 = np.pi/2)

    np.save("SkyImageData",ImageGrids)

def Minecraft():
    
    SkyData = np.load("himmelkule.npy")

    InputImage = Image.open("sample0000.png")
    FOV_PHI = np.deg2rad(90)
    FOV_THETA = np.deg2rad(90)

    PixelHeight, PixelWidth = 480,480

    # Generating Ranges
    RangeX = np.array(MaxMin(FOV_PHI))
    RangeY = -np.array(MaxMin(FOV_THETA))

    CoordinateGrid = GenerateCoordinateGrid(RangeX,RangeY,PixelWidth,PixelHeight)

    Angles = [
        (45,90),
        (135,90),
        (225,90),
        (315,90),
        (45,0),
        (135,0),
        (45,180),
        (135,180)
    ]

    for i in range(len(Angles)):

        phi = np.deg2rad(Angles[i][0])
        theta = np.deg2rad(Angles[i][1])

        AngleGrid = GenerateAngleGrid(CoordinateGrid,phi0 = phi,theta0 = theta,Wrap = False)

        ImageGrid = GenerateImage(SkyData,AngleGrid)
        # Reproducing sample2000.png
        OutputImage = Image.fromarray(ImageGrid)
        OutputImage.save(f"panorama_{i}.png") # Make new png

if __name__ == "__main__":

    if(len(sys.argv) == 1):
        print("No arguments. Running main...")
        Main()
    elif(sys.argv[1].lower() == "minecraft"):
        print("Minecraft time")
        Minecraft()
    else:
        print("No matching arguments. Running main...")
        Main()
    