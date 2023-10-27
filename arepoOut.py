# Importing libraries
import numpy as np

# Function to export the created particle data to a usable arepo file
def arepoOut(ngas):
    # Starting with the header file
    nPartArray = np.zeros(6, dtype=np.int32)

    # Setting the number of particles
    nPartArray[0] = ngas
    
