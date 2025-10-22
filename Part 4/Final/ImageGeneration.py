# Ikke brukt kodemal!!
# Skrevet av Bastian Eggum Huuse

import sys
import numpy as np
from numba import njit
from PIL import Image

import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

# This file contains many small functions that do small things
# This makes the structure of the program difficult to follow from the code alone
# Refer to the flowchart if this becomes the case.

@njit
def GetImageDimentions(Pixels):

    """
    Function that gets the dimentions of a given image

    Parameters:
    Pixels : Array(float) | Image to find dimentions of, converted to array

    returns:
    int | Height of image in pixels
    int | Width of image in pixels
    """
    
    Height = len(Pixels)
    Width = len(Pixels[0])

    return Height,Width

@njit
def MaxFunction(Alpha):

    """
    Function that finds maximum point along an axis, by using FOV along said axis

    Parameters:
    Alpha : float | FOV along the axis we want to find max along

    returns:
    float | max point along the same axis as the FOV
    """
    
    Value = (2*np.sin(Alpha/2))/(1+np.cos(Alpha/2))
    return Value

@njit
def MaxMin(FOV):

    """
    Function that finds maximum and minimum points along an axis, by using FOV along said axis

    Parameters:
    Alpha : Array(float) |  FOV along the axis we want to find max and min along

    returns:
    float | min point along the same axis as the FOV
    float | max point along the same axis as the FOV
    """

    Value = MaxFunction(FOV)

    Max = Value
    Min = -Value

    return Min,Max

@njit
def Rho(x,y):

    """Support function that describes the property Rho in the formulas for Phi and Theta"""

    return (x**2 + y**2)**(1/2)

@njit
def Beta(x,y):
    
    """Support function that describes the property Beta in the formulas for Phi and Theta"""

    rho = Rho(x,y)
    return 2* np.arctan(rho/2)

@njit
def Phi(phi0,theta0,x,y):

    """
    Formula that finds corresponding Phi value for a given x and y value

    Parameters:
    phi0   : float | Base phi in radians
    theta0 : float | Base theta in radians
    x      : float | x position along sphere
    y      : float | y position along sphere

    returns : float | Phi corresponding to given (x,y) position
    """

    rho = Rho(x,y)
    beta = Beta(x,y)

    return phi0 + np.arctan((x * np.sin(beta))/(rho*np.sin(theta0) * np.cos(beta) - y * np.cos(theta0)*np.sin(beta)))

@njit
def Theta(phi0,theta0,x,y):

    """
    Formula that finds corresponding Theta value for a given x and y value

    Parameters:
    phi0   : float | Base phi in radians
    theta0 : float | Base theta in radians
    x      : float | x position along sphere
    y      : float | y position along sphere

    returns : float | Theta corresponding to given (x,y) position
    """

    rho = Rho(x,y)
    beta = Beta(x,y)

    return theta0 - np.arcsin(np.cos(beta)*np.cos(theta0) + (y/rho)*np.sin(beta)*np.sin(theta0))

@njit
def GenerateCoordinateGrid(RangeX,RangeY,PixelWidth,PixelHeight):

    """
    Function that generates a grid of cartesian coordinates.
    Each pixel on the sample image has a corresponding coordinate.

    Parameters:
    RangeX      : tuple | Tuple containing min and max X values
    RangeY      : tuple | Tuple containing min and max Y values
    PixelWidth  : int   | Width of sample image in pixels
    PixelHeight : int   | Height of sample image in pixels

    returns : Array(float) | Grid of cartesian coordinates
    """

    # Generating linearly spaced values between the minimums and maximums
    ValuesX = np.linspace(RangeX[0],RangeX[1],PixelWidth)
    ValuesY = np.linspace(RangeY[0],RangeY[1],PixelHeight)

    # The image gives us height as the first dimention, but we want width as the first dimention
    CoordinateGrid = np.zeros((PixelWidth,PixelHeight,2))

    # Looping through grid, and assigning corresponding values
    for ix in range(len(ValuesX)):
        for iy in range(len(ValuesY)):
            CoordinateGrid[ix][iy] = np.array([ValuesX[ix],ValuesY[iy]])

    return CoordinateGrid

@njit
def GenerateAngleGrid(CoordinateGrid,phi0,theta0,Wrap = False):
    
    """
    Function that generates a grid of Spherical coordinates (without the radius),
    where each pixel of the sample image has one coordinate.

    Parameters:
    CoordinateGrid : Array(float) | Grid of cartesian coordinates
    phi0   : float | Base phi in radians
    theta0 : float | Base theta in radians
    Wrap   : bool  | Wether the function should constrict phi and theta (unused)
    """

    # Creating an empty grid
    AngleGrid = np.zeros((len(CoordinateGrid),len(CoordinateGrid[0]),2))

    # Looping through grid
    for ix in range(len(CoordinateGrid)):
        for iy in range(len(CoordinateGrid[ix])):
            
            # Finding x and y for this point in the grid
            x = CoordinateGrid[ix][iy][0]
            y = CoordinateGrid[ix][iy][1]

            # Computing theta and phi using x and y
            theta = Theta(phi0,theta0,x,y)
            phi = Phi(phi0,theta0,x,y)

            if(Wrap):
                if theta > np.pi:
                    theta = np.pi
                while theta < 0:
                    theta = 0

            # Assigning this point in the grid a theta,phi pair.
            AngleGrid[ix][iy] = np.array([theta,phi])

    return AngleGrid

def GenerateImage(SkyData,AngleGrid):
    
    """
    Method that generates an array that can be turned into an image.
    This is done through the method get_sky_image_pixel from SpaceMission

    Parameters:
    SkyData   : Array(float) | Array Containing the RGB values for all possible points on the sky sphere
    AngleGrid : Array(float) | grid of Theta,phi pairs, each point on the grid corresponding to a pixel.

    returns : Array(float) | Array Containing the RGB values of all pixels corresponding to our Angle Grid
    """

    # Notice the third dimention is 3 this time, since we want RGB
    ImageGrid = np.zeros((AngleGrid.shape[1],AngleGrid.shape[0],3))

    # Looping through all points on the grid
    for ix in range(len(AngleGrid)):
        for iy in range(len(AngleGrid[0])):
            # Getting the index of the current pixel (using the angle grid)
            PixelIndex = SpaceMission.get_sky_image_pixel(AngleGrid[ix][iy][0],AngleGrid[ix][iy][1])
            # Getting the current RGB values (The first two elements in each Skydata point are unused)
            # Also note that since images want the height first, we also flip the indexing here
            ImageGrid[iy][ix] = SkyData[PixelIndex][2:]

    # Returning the image with the correct datatyping
    return ImageGrid.astype('uint8')

def GenerateImageRange(Skydata,CoordinateGrid,theta0):

    """
    Generating 360 images (one for each degree), that are snapshots of one entire rotation of phi.

    Parameters:
    SkyData        : Array(float) | Array Containing the RGB values for all possible points on the sky sphere
    CoordinateGrid : Array(float) | grid of x,y pairs, each point on the grid corresponding to a pixel. (This is the same for every image)
    theta0         : float        | Base theta angle in radians (This is the same for every image)
    
    returns : Array(float) | Array Containing all 360 images in Array form.
    """

    # Creating an array of all angles
    phi0 = np.arange(0,360,1).astype("float32")

    # Creating a list to contain all 360 images
    Images = np.zeros((len(phi0),CoordinateGrid.shape[1],CoordinateGrid.shape[0],3),dtype="uint8")

    # Looping over every angle
    for i in range(len(phi0)):

        phi0[i] = np.deg2rad(phi0[i])

        # Generating and saving one image
        AngleGrid = GenerateAngleGrid(CoordinateGrid,phi0[i],theta0)
        ImageGrid = GenerateImage(Skydata,AngleGrid)
        Images[i] = ImageGrid

    return Images

def Main():

    # Loading the file containing all sky pixels
    SkyData = np.load("himmelkule.npy")

    # Loading out sample image, and defining our FOV
    InputImage = Image.open("sample0000.png")
    FOV_PHI = np.deg2rad(70)
    FOV_THETA = np.deg2rad(70)

    # Processing image (Turning it into an array)
    InputImage = InputImage
    Pixels = np.array(InputImage)
    PixelHeight, PixelWidth = GetImageDimentions(Pixels)

    # Generating Ranges (note that RangeY is reversed, since Theta goes down when Y goes up!)
    RangeX = np.array(MaxMin(FOV_PHI))
    RangeY = -np.array(MaxMin(FOV_THETA))

    # Generating grid of x and y values:
    CoordinateGrid = GenerateCoordinateGrid(RangeX,RangeY,PixelWidth,PixelHeight)
    # Generating grid of theta and phi values:
    AngleGrid = GenerateAngleGrid(CoordinateGrid,phi0 = 0,theta0 = np.pi/2)
    # Generating image based on these angles
    ImageGrid = GenerateImage(SkyData,AngleGrid)

    # Reproducing sample2000.png
    OutputImage = Image.fromarray(ImageGrid)
    OutputImage.save("Output0000.png") # Make new png

    # running for all phi0 between 0 and 2pi
    ImageGrids = GenerateImageRange(SkyData,CoordinateGrid,theta0 = np.pi/2)

    # Saving all generated images to files
    np.save("SkyImageData",ImageGrids)

def Minecraft():

    """
    This method is not relevant to the task, and was just done for fun. Therefore, it isn't properly documented :)
    """

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

    # Under normal circumstances, we simply run the function Main()
    if(len(sys.argv) == 1):
        print("No arguments. Running main...")
        Main()
        
    elif(sys.argv[1].lower() == "minecraft"):
        print("Minecraft time")
        Minecraft()
    else:
        print("No matching arguments. Running main...")
        Main()

"""
String that runs code: python ImageGeneration.py

Output:

No arguments. Running main...
"""
    