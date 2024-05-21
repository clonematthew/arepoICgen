# Importing libraries
import numpy as np

# Constants
gasConstant = 8.31

# Defining mass of the particles 
def masses(ngas, totalMass, verbose=False):
    # Calculate mass of each particle 
    particleMass = totalMass / ngas

    # Set array of particle masses
    pMass = np.ones(ngas, dtype=np.float64) * particleMass

    if verbose:
        # Printing the mass
        print("Total desired mass: %s" % totalMass)
        print("Initial particle mass: {:.5f}".format(particleMass))

    return pMass

# Defining energy of the particles
def thermalEnergy(ngas, temperature, mu, verbose=False):
    # Calculating internal energy
    energy = (3./2.) * temperature * gasConstant / mu 

    # Calculating a sound speed
    cs = np.sqrt(energy * 2./3.)

    # Allocating this energy to each particle
    pEnergy = np.ones(ngas, dtype=np.float64) * energy

    if verbose:
        # Printing values
        print("Initial particle energy: {:.2f}".format(energy))
        print("Initial sound speed: {:.2f}".format(cs))

    return pEnergy
