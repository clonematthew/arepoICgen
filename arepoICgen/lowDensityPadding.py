# Importing libraries
import numpy as np
from random import random

# Function to pad box with low density particles for a spherical cloud
def padBox(ngas, pos, vels, pMass, pIDs, pEnergy, boxDims, tempFactor):
    # Use 2% of the number of particles to pad the box 
    nPaddingParticles = int(0.02 * ngas)

    print("Padding the box with %s new particles" % nPaddingParticles)

    # Create new arrays that are long enough for all the particles
    newPos = np.zeros((3, nPaddingParticles), dtype=np.float64)
    newVels = np.zeros((3, nPaddingParticles), dtype=np.float64)
    newMass = np.zeros(nPaddingParticles, dtype=np.float64)
    newIDs = np.zeros(nPaddingParticles, dtype=np.int32)
    newEnergy = np.zeros(nPaddingParticles, dtype=np.float64)
    newRho = np.zeros(nPaddingParticles, dtype=np.float64)

    # Storing the old dimensions of our cloud
    xminOld = np.min(pos[0])
    xmaxOld = np.max(pos[0])
    yminOld = np.min(pos[1])
    ymaxOld = np.max(pos[1])
    zminOld = np.min(pos[2])
    zmaxOld = np.max(pos[2])

    # Append these new arrays to the end of our old ones
    pos = np.append(pos, newPos, axis=1)
    vels = np.append(vels, newVels, axis=1)
    pMass = np.append(pMass, newMass)
    pIDs = np.append(pIDs, newIDs)
    pEnergy = np.append(pEnergy, newEnergy)
    pRho = np.append(np.ones(ngas), newRho)

    # Getting the centre point of each dimension
    xMid = 0.5 * (xmaxOld + xminOld)
    yMid = 0.5 * (ymaxOld + yminOld)
    zMid = 0.5 * (zmaxOld + zminOld)

    # Getting the limits of our box
    xmin = xMid - boxDims[0] * (xmaxOld + xminOld) / 2
    xmax = xMid + boxDims[0] * (xmaxOld + xminOld) / 2
    ymin = yMid - boxDims[1] * (ymaxOld + yminOld) / 2
    ymax = yMid + boxDims[1] * (ymaxOld + yminOld) / 2
    zmin = zMid - boxDims[2] * (zmaxOld + zminOld) / 2
    zmax = zMid + boxDims[2] * (zmaxOld + zminOld) / 2

    # Working out volumes
    boxVolume = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
    cloudVolume = (xmaxOld - xminOld) * (ymaxOld - yminOld) * (zmaxOld - zminOld)

    # Calculating the density of the cloud
    cloudMass = np.sum(pMass)
    cloudDensity = cloudMass / cloudVolume

    # Setting the density of the particles within the cloud to this
    pRho = pRho * cloudDensity

    # Calculating mass of the particles we'll pad with
    newParticleMass = (0.01 * cloudDensity) * (boxVolume - cloudVolume) / nPaddingParticles

    # Randomly spraying the particles around the box
    placedPoints = 0 

    while placedPoints < nPaddingParticles:
        # Trying an x, y and z point
        xTry = xmin + (xmax - xmin) * random()
        yTry = ymin + (ymax - ymin) * random()
        zTry = zmin + (zmax - zmin) * random()

        # Checking if they're outside our dimensions
        i = 0
        if xTry > xminOld and xTry < xmaxOld:
            i += 1
        if yTry > yminOld and yTry < ymaxOld:
            i += 1
        if zTry > zminOld and zTry < zmaxOld:
            i += 1

        if i == 3:
            pass
        else:
            placedPoints += 1
            pID = ngas + placedPoints -1

            # Placing the particles inside the arrays
            pos[0,pID] = xTry
            pos[1,pID] = yTry
            pos[2,pID] = zTry
            pIDs[pID] = pID + 1
            pEnergy[pID] = pEnergy[0] * tempFactor
            pMass[pID] = newParticleMass
            pRho[pID] = 0.01 * cloudDensity

    return pos, vels, pMass, pIDs, pEnergy, pRho, (ngas+nPaddingParticles)

# Function to pad a spherical grid
def padSphere(ngas, pos, vels, pMass, pIDs, pEnergy, boxDims, tempFactor):
    # Use 2% of the number of particles to pad the box 
    nPaddingParticles = int(0.02 * ngas)

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

    # Working out the centre of the cloud and the box size
    mtot = np.sum(pMass)
    xcom = np.sum(pos[0] * pMass) / mtot
    ycom = np.sum(pos[1] * pMass) / mtot
    zcom = np.sum(pos[2] * pMass) / mtot

    # Getting dimensions from this
    xmin = boxDims[0] * np.min(pos[0])
    xmax = boxDims[0] * np.max(pos[0])
    ymin = boxDims[1] * np.min(pos[1])
    ymax = boxDims[1] * np.max(pos[1])
    zmin = boxDims[2] * np.min(pos[2])
    zmax = boxDims[2] * np.max(pos[2])

    # Find radius of the cloud
    cloudRadius = np.max(np.sqrt((pos[0] - xcom)**2 + (pos[1] - ycom)**2 + (pos[2] - zcom)**2))
    cloudVolume = np.pi * cloudRadius**3 * 4./3.

    # Determining the density of the cloud
    cloudDensity = 3.*mtot / (4 * np.pi) / cloudRadius**3

    # Setting this as denisty of particles in the cloud
    pRho = pRho * cloudDensity

    # Calculating mass of the particles we'll pad with
    newParticleMass = (0.01 * cloudDensity) * cloudVolume / nPaddingParticles

    # Randomly spraying particles around the sphere
    placedPoints = 0

    while placedPoints < nPaddingParticles:
        # Trying an x, y and z point
        xTry = xmin + (xmax - xmin) * random()
        yTry = ymin + (ymax - ymin) * random()
        zTry = zmin + (zmax - zmin) * random()

        # Finding distance from the centre
        r = np.sqrt((xTry - xcom)**2 + (yTry - ycom)**2 + (zTry - zcom)**2)

        # Adding if its outside the cloud
        if r > cloudRadius:
            placedPoints += 1
            pID = ngas + placedPoints -1

            # Placing the particles inside the arrays
            pos[0,pID] = xTry
            pos[1,pID] = yTry
            pos[2,pID] = zTry
            pIDs[pID] = pID + 1
            pEnergy[pID] = pEnergy[0] * tempFactor
            pMass[pID] = newParticleMass
            pRho[pID] = 0.01 * cloudDensity
    
        else:
            pass

    return pos, vels, pMass, pIDs, pEnergy, pRho, (ngas+nPaddingParticles)

def padEllipse(ngas, pos, vels, pMass, pIDs, pEnergy, boxDims, eX, eY, eZ, tempFactor):
    # Convert the ellipse dimensions
    eX = 3.09e18 * eX
    eY = 3.09e18 * eY
    eZ = 3.09e18 * eZ

    # Use 2% of the number of particles to pad the box 
    nPaddingParticles = int(0.02 * ngas)

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

    # Working out the centre of the cloud and the box size
    mtot = np.sum(pMass)
    xcom = np.sum(pos[0] * pMass) / mtot
    ycom = np.sum(pos[1] * pMass) / mtot
    zcom = np.sum(pos[2] * pMass) / mtot

    # Getting dimensions from this
    maxDimension = np.max(pos)
    minDimensionX = np.min(pos) * boxDims[0]
    minDimensionY = np.min(pos) * boxDims[1]
    minDimensionZ = np.min(pos) * boxDims[2]
    maxDimensionX = boxDims[0] * maxDimension
    maxDimensionY = boxDims[1] * maxDimension
    maxDimensionZ = boxDims[2] * maxDimension
    boxVolume = (maxDimensionX - minDimensionX) * (maxDimensionY - minDimensionY) * (maxDimensionZ - minDimensionZ)

    # Getting the volume of the cloud
    cloudVolume = (np.pi * 4/3) * (eX * eY * eZ)
    cloudMass = np.sum(pMass)
    cloudDensity = cloudMass / cloudVolume

    # Setting the density of the particles within the cloud to this
    pRho = pRho * cloudDensity

    # Calculating mass of the particles we'll pad with
    newParticleMass = (0.01 * cloudDensity) * (cloudVolume) / nPaddingParticles

    # Randomly spraying the particles around the box
    placedPoints = 0 

    while placedPoints < nPaddingParticles:
        # Trying an x, y and z point
        xTry = minDimensionX + (maxDimensionX - minDimensionX) * random()
        yTry = minDimensionY + (maxDimensionY - minDimensionY) * random()
        zTry = minDimensionZ + (maxDimensionZ - minDimensionZ) * random()

        # Adding if its outside the cloud
        if abs(xTry) < eX and abs(yTry) < eY and abs(zTry) < eZ:
            pass
        else:
            placedPoints += 1
            pID = ngas + placedPoints -1

            # Placing the particles inside the arrays
            pos[0,pID] = xTry
            pos[1,pID] = yTry
            pos[2,pID] = zTry
            pIDs[pID] = pID + 1
            pEnergy[pID] = pEnergy[0] * tempFactor
            pMass[pID] = newParticleMass
            pRho[pID] = 0.01 * cloudDensity

    return pos, vels, pMass, pIDs, pEnergy, pRho, (ngas+nPaddingParticles)   

def padCylinder(ngas, pos, vels, pMass, pIDs, pEnergy, boxDims, length, radius, tempFactor):
    # Use 2% of the number of particles to pad the box 
    nPaddingParticles = int(0.02 * ngas)

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

    # Find the length of the cloud in each dimension
    xWidth = np.max(pos[0]) - np.min(pos[0])
    yWidth = np.max(pos[1]) - np.min(pos[1])
    zWidth = np.max(pos[2]) - np.min(pos[2])

    # Working out the final size of the cloud
    minDimension = np.min(pos)
    maxDimension = np.max(pos)

    minDimension *= boxDims
    maxDimension *= boxDims

    boxVolume = (maxDimension - minDimension)**3

    # Working out the density 
    cloudVolume = (np.pi * 4/3) * (xWidth * yWidth * zWidth)
    cloudDensity = cloudVolume / np.sum(pMass)

    # Assigning the density of the padding particles
    newParticleMass = (0.01 * cloudDensity) * (cloudVolume) / nPaddingParticles

    # Finding particles to pad the cloud with
    placedPoints = 0 

    while placedPoints < nPaddingParticles:
        # Trying an x, y and z point
        xTry = minDimension + (maxDimension - minDimension) * random()
        yTry = minDimension + (maxDimension - minDimension) * random()
        zTry = minDimension + (maxDimension - minDimension) * random()

        # Adding if its outside the cloud
        if abs(xTry) <= 0.5 * length * 3.09e18 and np.sqrt(yTry**2 + zTry**2) <= radius * 3.09e18:
            pass
        else:
            placedPoints += 1
            pID = ngas + placedPoints -1

            # Placing the particles inside the arrays
            pos[0,pID] = xTry
            pos[1,pID] = yTry
            pos[2,pID] = zTry
            pIDs[pID] = pID + 1
            pEnergy[pID] = pEnergy[0] * tempFactor
            pMass[pID] = newParticleMass
            pRho[pID] = 0.01 * cloudDensity

    return pos, vels, pMass, pIDs, pEnergy, pRho, (ngas+nPaddingParticles)  