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
    rBonnorEbert = 3.09e18 * 0.0043016 * (beMass/1.991e33) * 6.5 / ((cs/1e5)**2 * 15.85)
    print("Bonnor Ebert Radius: {:.2e}".format(rBonnorEbert))
     
    outsideRho = 1.66e-22
    centralDensity = 14.305 * outsideRho
    rCharacteristic = cs / (np.sqrt(4 * np.pi * G * centralDensity))

    #centralDensity = (6.5**2 / np.sqrt(4 * np.pi)) * (cs**2 / G) * (1/rBonnorEbert**2)
    print("Central Density: {:.2e}".format(centralDensity))
    #rCharacteristic = cs / np.sqrt(4 * np.pi * G * centralDensity)
    print("Characteristic Radius: {:.2e}".format(rCharacteristic))
    
    # Create bins of radius and corresponding density
    radiusBins = 10**np.linspace(0, np.log10(rBonnorEbert), 1000)
    densityBins = centralDensity * rCharacteristic**2 / (rCharacteristic**2 + radiusBins**2)

    # Find the centre of mass
    totMass = np.sum(mass)
    xC = np.sum(pos[0] * mass) / totMass
    yC = np.sum(pos[1] * mass) / totMass
    zC = np.sum(pos[2] * mass) / totMass
    
    # Find each particle's distance to the COM
    rCentre = np.sqrt((pos[0])**2 + (pos[1])**2 + (pos[2])**2)
    
    # Scale everything to the correct radius
    pos[0] = rBonnorEbert * pos[0] / np.max(rCentre) 
    pos[1] = rBonnorEbert * pos[1] / np.max(rCentre)
    pos[2] = rBonnorEbert * pos[2] / np.max(rCentre)
    rCentre = rBonnorEbert * rCentre / np.max(rCentre) 
    
    # Scale the mass of each particle
    for i in range(len(radiusBins)-1): 
        inRadius = np.where((rCentre > radiusBins[i]) & (rCentre < radiusBins[i+1]))
            
        # Work out how much mass should be in each shell
        shellMass = (4 * np.pi / 3) * (radiusBins[i+1]**3 - radiusBins[i]**3) * densityBins[i] 
        
        # Find the total mass of particles in this shell and scale
        cellMass = shellMass / len(inRadius[0])
        
        # Assign value
        mass[inRadius] = cellMass
        
    # Scale masses
    totMass = np.sum(mass)
    mass = mass * (beMass / totMass)
        
    # Work out cloud density
    volume = ((4 * np.pi / 3) * rBonnorEbert**3)
    cloudDensity = beMass / volume
    densityFrac = np.min(densityBins)/cloudDensity

    return mass, pos, densityFrac/15, volume/(3.09e18**3)