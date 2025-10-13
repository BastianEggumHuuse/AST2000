# Ikke brukt kodemal!!
# Skrevet av Bastian Eggum Huuse

import numpy as np

from numba import njit

from PIL import Image

import ast2000tools.constants as const
import ast2000tools.utils     as utils
from ast2000tools.space_mission import SpaceMission

if __name__ == "__main__":

    SkyImageData = np.load("SkyImageData.npy")

    i = 0
    for img in SkyImageData:
        OutputImage = Image.fromarray(img)
        OutputImage.save(f"grid{i}.png") # Make new png
        i += 1