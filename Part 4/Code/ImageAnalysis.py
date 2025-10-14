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

    Img1 = Image1.copy().astype("float32")
    Img2 = Image2.copy().astype("float32")

    Diff = 0

    # Comparing iteratively
    for i in range(len(Img1)):
        for j in range(len(Img1[i])):

            dR = Img1[i][j][0] - Img2[i][j][0]
            dG = Img1[i][j][1] - Img2[i][j][1]
            dB = Img1[i][j][2] - Img2[i][j][2]

            dC = (dR**2 + dG**2 + dB**2)
            Diff += dC

    return Diff

@njit
def CompareImageRange(InputImage,SkyImageData,N = 360):

    ComparisonImages = SkyImageData[:N]
    Diffs = np.zeros(N)

    # Calculating all differences
    for i in range(len(Diffs)):
        Diffs[i] = CompareImages(InputImage,ComparisonImages[i])

    # Finding the lowest diff (index)
    MinIndex = 0
    for i in range(1,len(Diffs)):

        if(Diffs[i] < Diffs[MinIndex]):
            MinIndex = i

    # Since we have 360 different angles indexed,
    # we can return the angle by simply returning the index
    # (Note that this gives angle in degrees)

    return MinIndex


if __name__ == "__main__":

    SkyImageData = np.load("SkyImageData.npy")

    InputPath = "sample0200.png"
    if(len(sys.argv) > 1):
        InputPath = sys.argv[1]

    print(f"Loading image {InputPath}")

    InputImage = Image.open(InputPath)
    InputImageArray = np.array(InputImage)

    print(f"Finding closest match through least squares...\n")

    phi = CompareImageRange(InputImageArray,SkyImageData,N = 359)

    print(f"Closest angle was {phi} degrees.")
    print(f"Showing input image along with image at {phi} degrees:")

    OutputImageArray = SkyImageData[phi]
    OutputImage = Image.fromarray(OutputImageArray)

    # Combining the two images
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