# Importing libraries
import numpy as np
from random import random

# Function to pad the box around the cloud
def padBox(ngas, pos, vels, pMass, pEnergy, volume, config, params):
    # Check if we are to pad the box
    if config["padBox"]:
        # Check if we are using the special case of a BE sphere
        if config["extras"] == "bonnorEbert":
            # We want the outer pressure to match the BE sphere
            innerPressure = np.min(pMass[pMass > 0]) * 8.31 * params["temp"] / (params["mu"])
            outerDensity = innerPressure * params["mu"] / (8.31 * params["tempFactor"] * params["temp"])
            params["paddingDensity"] = outerDensity / np.min(pMass[pMass > 0])

        # Pad the box 
        ngas, pos, vels, pMass, pEnergy = padGeneric(ngas, pos, vels, pMass, pEnergy, volume, config, params)            
        
        return ngas, pos, vels, pMass, pEnergy  
            
    else:
        return ngas, pos, vels, pMass, pEnergy 
                    
# Generic function for the padding
def padGeneric(ngas, pos, vels, pMass, pEnergy, volume, config, params):
    # Use 2% of the number of particles to pad the box 
    nPaddingParticles = int(params["paddingPercent"] * ngas)

    if config["verbose"]:
        print("Padding the box with %s new particles" % nPaddingParticles)

    # Create new arrays that are long enough for all the particles
    newPos = np.zeros((3, nPaddingParticles), dtype=np.float64)
    newVels = np.zeros((3, nPaddingParticles), dtype=np.float64)
    newMass = np.zeros(nPaddingParticles, dtype=np.float64)
    newEnergy = np.zeros(nPaddingParticles, dtype=np.float64)

    # Append these new arrays to the end of our old ones
    pos = np.append(pos, newPos, axis=1)
    vels = np.append(vels, newVels, axis=1)
    pMass = np.append(pMass, newMass)
    pEnergy = np.append(pEnergy, newEnergy)

    # Find the bounds of the cloud
    xmin = np.min(pos[0])
    xmax = np.max(pos[0])
    ymin = np.min(pos[1])
    ymax = np.max(pos[1])
    zmin = np.min(pos[2])
    zmax = np.max(pos[2])

    # Find the centre of mas of the cloud
    xcom = np.sum(pMass * pos[0]) / np.sum(pMass)
    ycom = np.sum(pMass * pos[1]) / np.sum(pMass)
    zcom = np.sum(pMass * pos[2]) / np.sum(pMass)

    # Getting dimensions of the new box
    maxDimension = np.max(pos)
    minDimension = np.min(pos)
    minDimensionX = params["boxSize"][0] * minDimension
    minDimensionY = params["boxSize"][1] * minDimension
    minDimensionZ = params["boxSize"][2] * minDimension
    maxDimensionX = params["boxSize"][0] * maxDimension
    maxDimensionY = params["boxSize"][1] * maxDimension
    maxDimensionZ = params["boxSize"][2] * maxDimension

    # Working out volume of the box
    boxVolume = (maxDimensionX - minDimensionX) * (maxDimensionY - minDimensionY) * (maxDimensionZ - minDimensionZ)

    # Work out the cloud volume and denisty 
    cloudVolume = volume #* (3.09e18)**3
    cloudMass = np.sum(pMass)
    cloudDensity = cloudMass / cloudVolume
    cloudDimensions = [xmin, xmax, ymin, ymax, zmin, zmax]
    cloudCentre = [xcom, ycom, zcom]
    cloudRadius = np.max([(xmax-xmin)/2, (ymax-ymin)/2, (zmax-zmin)/2])

    pc = 3.09e18
    if config["verbose"]:
        print("Cloud Density of {:.2e}".format(cloudDensity))
        print("Cloud Centre at {:.2f},{:.2f},{:.2f}".format(xcom/pc, ycom/pc, zcom/pc))
        print("Cloud has Dimensions: {:.2f}-{:.2f},{:.2f}-{:.2f},{:.2f}-{:.2f}".format(cloudDimensions[0]/pc, cloudDimensions[1]/pc, cloudDimensions[2]/pc, cloudDimensions[3]/pc, cloudDimensions[4]/pc, cloudDimensions[5]/pc))
        print("Box has Dimensions: {:.2f}-{:.2f},{:.2f}-{:.2f},{:.2f}-{:.2f}".format(minDimensionX/pc, maxDimensionX/pc, minDimensionY/pc, maxDimensionY/pc, minDimensionZ/pc, maxDimensionZ/pc))

    # Calculating mass of the particles we'll pad with 
    newParticleMass = (params["paddingDensity"] * cloudDensity) * (boxVolume-cloudVolume) / nPaddingParticles 
    minCellDensity = np.min(pMass[pMass > 0])

    # Randomly spraying the particles around the box
    placedPoints = 0 

    while placedPoints < nPaddingParticles:
        # Trying an x, y and z point
        xTry = minDimensionX + (maxDimensionX - minDimensionX) * random()
        yTry = minDimensionY + (maxDimensionY - minDimensionY) * random()
        zTry = minDimensionZ + (maxDimensionZ - minDimensionZ) * random()

        # Checking if they're outside our dimensions
        if outsideBox(config["grid"], xTry, yTry, zTry, cloudDimensions, cloudCentre, cloudRadius):
            placedPoints += 1
            pID = ngas + placedPoints -1

            # Placing the particles inside the arrays
            pos[0,pID] = xTry
            pos[1,pID] = yTry
            pos[2,pID] = zTry
            pEnergy[pID] = pEnergy[0] * params["tempFactor"]
            
            if config["extras"] == "bonnorEbert":
                pMass[pID] = minCellDensity * params["paddingDensity"]
            else:
                pMass[pID] = newParticleMass 
    
    return  (ngas+nPaddingParticles), pos, vels, pMass, pEnergy

# Check if the particle is inside each shape's geometry
def outsideBox(gridType, xTry, yTry, zTry, cloudDimensions, cloudCentre, cloudRadius):
    # Box geometry setups
    if gridType == "boxGrid" or gridType == "boxRan":
        i = 0
        if xTry > cloudDimensions[0] and xTry < cloudDimensions[1]:
            i += 1
        if yTry > cloudDimensions[2] and yTry < cloudDimensions[3]:
            i += 1
        if zTry > cloudDimensions[4] and zTry < cloudDimensions[5]:
            i += 1

        if i == 3:
            return False
        else:
            return True
        
    # Spherical geometry setups
    elif gridType == "sphereGrid" or gridType == "sphereRan":
        r = np.sqrt((xTry - cloudCentre[0])**2 + (yTry - cloudCentre[1])**2 + (zTry - cloudCentre[2])**2)

        if r <= cloudRadius:
            return False
        else:
            return True
        
    # Ellipsoidal geometry setups
    elif gridType == "ellipseRan":
        xx = ((cloudDimensions[1] - cloudDimensions[0])/2)**2 
        yy = ((cloudDimensions[3] - cloudDimensions[2])/2)**2 
        zz = ((cloudDimensions[5] - cloudDimensions[4])/2)**2 

        if (xTry*xTry/xx + yTry*yTry/yy + zTry*zTry/zz) <= 1:
            return False
        else:
            return True

    # Cylindrical geometry setups:
    elif gridType == "cylinderRan":
        length = cloudDimensions[1] - cloudDimensions[0]
        radius = (cloudDimensions[3] - cloudDimensions[2]) / 2

        if abs(xTry) <= 0.5 * length  and np.sqrt(yTry**2 + zTry**2) <= radius: 
            return False
        else:
            return True
        
    # Unimplemented setups
    else:
        print("Geometry not implemented!")
        return False        