# Importing libraries
import numpy as np
from random import random 

# Generic function for the padding
def padGeneric(ngas, pos, vels, pMass, pIDs, pEnergy, volume, boxDims, gridType, tempFactor, paddingPercent=0.02, padDensity=0.01, verbose=False):
    # Use 2% of the number of particles to pad the box 
    nPaddingParticles = int(paddingPercent * ngas)

    if verbose:
        print("Padding the box with %s new particles" % nPaddingParticles)

    # Create new arrays that are long enough for all the particles
    newPos = np.zeros((3, nPaddingParticles), dtype=np.float64)
    newVels = np.zeros((3, nPaddingParticles), dtype=np.float64)
    newMass = np.zeros(nPaddingParticles, dtype=np.float64)
    newIDs = np.zeros(nPaddingParticles, dtype=np.int32)
    newEnergy = np.zeros(nPaddingParticles, dtype=np.float64)
    newRho = np.zeros(nPaddingParticles, dtype=np.float64)

    # Append these new arrays to the end of our old ones
    pos = np.append(pos, newPos, axis=1)
    vels = np.append(vels, newVels, axis=1)
    pMass = np.append(pMass, newMass)
    pIDs = np.append(pIDs, newIDs)
    pEnergy = np.append(pEnergy, newEnergy)
    pRho = np.append(np.ones(ngas), newRho)

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
    minDimensionX = boxDims[0] * minDimension
    minDimensionY = boxDims[1] * minDimension
    minDimensionZ = boxDims[2] * minDimension
    maxDimensionX = boxDims[0] * maxDimension
    maxDimensionY = boxDims[1] * maxDimension
    maxDimensionZ = boxDims[2] * maxDimension

    # Working out volume of the box
    boxVolume = (maxDimensionX - minDimensionX) * (maxDimensionY - minDimensionY) * (maxDimensionZ - minDimensionZ)

    # Work out the cloud volume and denisty 
    cloudVolume = volume * (3.09e18)**3
    cloudMass = np.sum(pMass)
    cloudDensity = cloudMass / cloudVolume
    cloudDimensions = [xmin, xmax, ymin, ymax, zmin, zmax]
    cloudCentre = [xcom, ycom, zcom]
    cloudRadius = np.max([(xmax-xmin)/2, (ymax-ymin)/2, (zmax-zmin)/2])

    pc = 3.09e18
    if verbose:
        print("Cloud Density of {:.2e}".format(cloudDensity))
        print("Cloud Centre at {:.2f},{:.2f},{:.2f}".format(xcom/pc, ycom/pc, zcom/pc))
        print("Cloud has Dimensions: {:.2f}-{:.2f},{:.2f}-{:.2f},{:.2f}-{:.2f}".format(cloudDimensions[0]/pc, cloudDimensions[1]/pc, cloudDimensions[2]/pc, cloudDimensions[3]/pc, cloudDimensions[4]/pc, cloudDimensions[5]/pc))
        print("Box has Dimensions: {:.2f}-{:.2f},{:.2f}-{:.2f},{:.2f}-{:.2f}".format(minDimensionX/pc, maxDimensionX/pc, minDimensionY/pc, maxDimensionY/pc, minDimensionZ/pc, maxDimensionZ/pc))

    # Setting the density of the particles within the cloud to this
    pRho = pRho * cloudDensity

    # Calculating mass of the particles we'll pad with 
    newParticleMass = (padDensity * cloudDensity) * (boxVolume-cloudVolume) / nPaddingParticles 

    # Randomly spraying the particles around the box
    placedPoints = 0 

    while placedPoints < nPaddingParticles:
        # Trying an x, y and z point
        xTry = minDimensionX + (maxDimensionX - minDimensionX) * random()
        yTry = minDimensionY + (maxDimensionY - minDimensionY) * random()
        zTry = minDimensionZ + (maxDimensionZ - minDimensionZ) * random()

        # Checking if they're outside our dimensions
        if outsideBox(gridType, xTry, yTry, zTry, cloudDimensions, cloudCentre, cloudRadius):
            placedPoints += 1
            pID = ngas + placedPoints -1

            # Placing the particles inside the arrays
            pos[0,pID] = xTry
            pos[1,pID] = yTry
            pos[2,pID] = zTry
            pIDs[pID] = pID + 1
            pEnergy[pID] = pEnergy[0] * tempFactor
            pMass[pID] = newParticleMass 
            pRho[pID] = padDensity * cloudDensity

    return pos, vels, pMass, pIDs, pEnergy, pRho, (ngas+nPaddingParticles)  

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