# Importing libraries
import numpy as np
from random import random

# Function to pad box with low density particles for a spherical cloud
def padBox(ngas, pos, vels, pMass, pIDs, pEnergy, boxDims, tempFactor):
    # Use 2% of the number of particles to pad the box 
    nPaddingParticles = int(0.02 * ngas)

    # Create new arrays that are long enough for all the particles
    newPos = np.zeros((3, nPaddingParticles))
    newVels = np.zeros((3, nPaddingParticles))
    newMass = np.zeros(nPaddingParticles)
    newIDs = np.zeros(nPaddingParticles)
    newEnergy = np.zeros(nPaddingParticles)
    newRho = np.zeros(nPaddingParticles)

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

    # Setting new cloud dimensions
    xmin = - 0.5 * boxDims[0]
    xmax = 0.5 * boxDims[0]
    ymin = - 0.5 * boxDims[1]
    ymax = 0.5 * boxDims[1]
    zmin = - 0.5 * boxDims[2]
    zmax = 0.5* boxDims[2]

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
            pID = ngas + placedPoints - 1

            # Placing the particles inside the arrays
            pos[0,pID] = xTry
            pos[1,pID] = yTry
            pos[2,pID] = zTry
            pIDs[pID] = pID
            pEnergy[pID] = pEnergy[0] * tempFactor
            pMass[pID] = newParticleMass
            pRho[pID] = 0.01 * cloudDensity

    return pos, vels, pMass, pIDs, pEnergy, pRho

# Function to pad a spherical grid
def padSphere(ngas, pos, vels, pMass, pIDs, pEnergy, tempFactor):
    # Use 2% of the number of particles to pad the box 
    nPaddingParticles = int(0.02 * ngas)

    # Create new arrays that are long enough for all the particles
    newPos = np.zeros((3, nPaddingParticles))
    newVels = np.zeros((3, nPaddingParticles))
    newMass = np.zeros(nPaddingParticles)
    newIDs = np.zeros(nPaddingParticles)
    newEnergy = np.zeros(nPaddingParticles)
    newRho = np.zeros(nPaddingParticles)

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
    xmin = 2 * np.min(pos[0])
    xmax = 2 * np.max(pos[0])
    ymin = 2 *np.min(pos[1])
    ymax = 2* np.max(pos[1])
    zmin = 2 *np.min(pos[2])
    zmax = 2* np.max(pos[2])

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
            pID = ngas + placedPoints - 1

            # Placing the particles inside the arrays
            pos[0,pID] = xTry
            pos[1,pID] = yTry
            pos[2,pID] = zTry
            pIDs[pID] = pID
            pEnergy[pID] = pEnergy[0] * tempFactor
            pMass[pID] = newParticleMass
            pRho[pID] = 0.01 * cloudDensity
    
        else:
            pass

    return pos, vels, pMass, pIDs, pEnergy, pRho




