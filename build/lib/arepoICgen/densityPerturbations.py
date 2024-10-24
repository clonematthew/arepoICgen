# Importing libraries
import numpy as np

# Some useful constants
G = 6.67e-8
kB = 1.38e-16

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

    return pos, mass

def densityGradient(pos, mass, lowerDensity=0.66, upperDensity=1.33):
    distance = pos[0] + np.max(pos[0])
    return lowerDensity * mass + (upperDensity-lowerDensity) * mass * (distance/np.max(distance))

def bonnorEbert(ngas, pos, mass, temp, mu, beMass):
    # Calculate the sound speed
    cs = np.sqrt(kB * temp / (mu * 1.66e-24))
    beMass = beMass * 1.991e33
    
    # Calculate characteristic quantities 
    rBonnorEbert = G * beMass / (2.42 * cs**2)
    print("Bonnor Ebert Radius: {:.2e}".format(rBonnorEbert))
    centralDensity = (6.5**2 / np.sqrt(4 * np.pi)) * (cs**2 / G) * (1/rBonnorEbert**2)
    print("Central Density: {:.2e}".format(centralDensity))
    rCharacteristic = cs / np.sqrt(4 * np.pi * G * centralDensity)
    print("Characteristic Radius: {:.2e}".format(rCharacteristic))
    
    # Find the centre of mass
    totMass = np.sum(mass)
    xC = np.sum(pos[0] * mass) / totMass
    yC = np.sum(pos[1] * mass) / totMass
    zC = np.sum(pos[2] * mass) / totMass
    
    # Find each particle's distance to the COM
    rCentre = np.sqrt((pos[0] - xC)**2 + (pos[1] - yC)**2 + (pos[2] - zC)**2)
    
    # Scale everything to the correct radius
    pos[0] = pos[0] * (np.max(rCentre) / rBonnorEbert)
    pos[1] = pos[1] * (np.max(rCentre) / rBonnorEbert)
    pos[2] = pos[2] * (np.max(rCentre) / rBonnorEbert)
    rCentre = rCentre * (np.max(rCentre) / rBonnorEbert)
    
    # Apply the density profile
    mass = mass * rCharacteristic**2 / (rCharacteristic**2 + rCentre**2)
    
    # Work out the mass of the sphere
    totMass = np.sum(mass)
    
    # Adjust to fit the desired mass
    mass = mass * (beMass/totMass)
    
    return mass, pos