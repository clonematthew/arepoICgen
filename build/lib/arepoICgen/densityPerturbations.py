# Importing libraries
import numpy as np

# Some useful constants
G = 6.67e-8

# Function to apply a BB density fluctuation
def bossBodenheimer(ngas, pos, mass):
    # Calculate the centre of mass
    totMass = np.sum(mass)
    xCom = np.sum(pos[0] * mass) / totMass
    yCom = np.sum(pos[1] * mass) / totMass

    # Apply the density perturbation
    for i in range(ngas):
        # Find relative positions
        x = xCom - pos[0,i]
        y = yCom - pos[1,i]

        # Work out the angle 
        phi = np.arctan2(y, x)

        # Work out what the mass should be here
        mass[i] = mass[i] * (1 + 0.5 * np.cos(2*phi))
        
    print("Boss-Bodenheimer Density Perturbation Applied")

    return pos, mass

def densityGradient(pos, mass, lowerDensity=0.66, upperDensity=1.33):
    distance = pos[0] + np.max(pos[0])
    return lowerDensity * mass + (upperDensity-lowerDensity) * mass * (distance/np.max(distance))
