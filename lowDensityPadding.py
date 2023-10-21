# Importing libraries
import numpy as np

# Function to pad box with low density particles for a spherical cloud
def padBoxSphere(ngas, pos, vels, pMass, pIDs, pEnergy):
    # Use 2% of the number of particles to pad the box 
    nPaddingParticles = int(0.02 * ngas)
    ngasNew = ngas + nPaddingParticles

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

    # Working out the cloud centre and the box size
    






    return pos, vels, pMass, pIDs, pEnergy, pRho


