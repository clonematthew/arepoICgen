# Importing libraries
import numpy as np

# Function to apply a BB density fluctuation
def bossBodenheimer(ngas, pos, mass):
    # Calculate the centre of mass
    totMass = np.sum(mass)
    xCom = np.sum(pos[0] * mass) / totMass
    yCom = np.sum(pos[1] * mass) / totMass

    # Get the mass resolution
    mRes = totMass / ngas

    # Apply the density perturbation
    for i in range(ngas):
        # Find relative positions
        x = xCom - pos[0,i]
        y = yCom - pos[1,i]

        # Work out the angle 
        phi = np.arctan(y, x)

        # Work out what the mass should be here
        mass[i] = mRes * (1 + 0.5 * np.cos(2*phi))

    print("Boss-Bodenheimer Density Perturbation Applied")
    return pos, mass