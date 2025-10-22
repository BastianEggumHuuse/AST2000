# Ikke brukt kodemal!!
# Skrevet av Bastian Eggum Huuse

import sys
import numpy as np
from numba import njit

from PIL import Image,ImageFont, ImageDraw

import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

@njit
def CompareImages(Image1, Image2):

    '''
    Method that compares two images by taking the sum of the squares of the differences of their RGB values

    Parameters:
    Image1 : Array(float) | the first image of the two we want to compare, converted to an array
    Image2 : Array(float) | the second image of the two we want to compare, converted to an array

    returns : float       | sum of the squares of the differences of the RGB values of the two images
    '''

    # Copying arrays so we don't modify the originals
    Img1 = Image1.copy().astype("float32")
    Img2 = Image2.copy().astype("float32")

    Diff = 0

    # Comparing iteratively (Going through each pixel)
    for i in range(len(Img1)):
        for j in range(len(Img1[i])):

            # Computing the differences of the R, G, and B values
            dR = Img1[i][j][0] - Img2[i][j][0]
            dG = Img1[i][j][1] - Img2[i][j][1]
            dB = Img1[i][j][2] - Img2[i][j][2]

            # Summing the squares of these
            dC = (dR**2 + dG**2 + dB**2)
            
            # Adding this sum to the total difference
            Diff += dC

    return Diff

@njit
def CompareImageRange(InputImage,SkyImageData,N = 360):

    '''
    Comparing an input image with 360 other images, one for each degree.
    The method then finds the image with the lowest total difference, which is the closest match

    Parameters:
    InputImage   : Array(float) | the image we want to approximate, converted to an array
    SkyImageData : Array(float) | an array containing 360 images, one for each degree (0-359)
    N            : int          | how many of the degrees we want to compare InputImage to (0 <= N < 360)
    '''

    # Getting the N first images
    ComparisonImages = SkyImageData[:N]
    Diffs = np.zeros(N)

    # Calculating all differences
    for i in range(len(Diffs)):
        Diffs[i] = CompareImages(InputImage,ComparisonImages[i])

    # Finding the lowest difference (index)
    MinIndex = 0
    for i in range(1,len(Diffs)):

        if(Diffs[i] < Diffs[MinIndex]):
            MinIndex = i

    # Since we have 360 different angles indexed,
    # we can return the angle by simply returning the index
    # (Note that this gives angle in degrees)

    return MinIndex


if __name__ == "__main__":

    # Loading SkyImageData
    # This contains 360 images, one for each degree
    SkyImageData = np.load("SkyImageData.npy")

    # Collecting a sample image we want to find the angle of
    InputPath = "sample0435.png"
    if(len(sys.argv) > 1): # Path can be input as command line argument
        InputPath = sys.argv[1]

    print(f"Loading image {InputPath}")

    # Loading the image with PIL
    InputImage = Image.open(InputPath)
    InputImageArray = np.array(InputImage)

    print(f"Finding closest match through least squares...\n")

    # Comparing images
    phi = CompareImageRange(InputImageArray,SkyImageData,N = 359)

    print(f"Closest angle was {phi} degrees.")
    print(f"Showing input image along with image at {phi} degrees:")

    # Getting the array associated with the output angle phi, 
    # and turning this back into an image
    OutputImageArray = SkyImageData[phi]
    OutputImage = Image.fromarray(OutputImageArray)



    # The rest of this program is not part of the task, but is just a way to check if the task was
    # performed correctly. We combine the two images (input and output), and show them side by side.
    buffX = 20

    CombinedImage = Image.new("RGBA", (640 * 2 + buffX, 480),color=(255,255,255,255))
    CombinedImage.paste(InputImage,(0,0))
    CombinedImage.paste(OutputImage,(640 + buffX,0))

    # Adding text to image
    textBuff = 10
    font = ImageFont.truetype("arial.ttf", 28)
    ImageDraw.Draw(CombinedImage).text(
        (textBuff,0),
        "Input Image",
        (255,0,0),
        font = font
    )
    ImageDraw.Draw(CombinedImage).text(
        (640 + buffX + textBuff,0),
        "Output Image",
        (255,0,0),
        font = font
    )

    CombinedImage.show()

"""
String that runs code: python ImageAnalysis.py

This program also shows an image.

Output:
Loading image sample0435.png
Finding closest match through least squares...

Closest angle was 44 degrees.
Showing input image along with image at 44 degrees:
"""