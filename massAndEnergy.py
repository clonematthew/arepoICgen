# Importing libraries
import numpy as np

# Constants
gasConstant = 8.31

# Defining mass of the particles 
def masses(ngas, totalMass):
    # Calculate mass of each particle 
    particleMass = totalMass / ngas

    # Set array of particle masses
    pMass = np.ones(ngas) * particleMass

    return pMass

# Defining energy of the particles
def thermalEnergy(ngas, temperature, mu):
    # Calculating internal energy
    energy = (3./2.) * temperature * gasConstant / mu 

    # Calculating a sound speed
    cs = np.sqrt(energy * 2./3.)

    # Allocating this energy to each particle
    pEnergy = np.ones(ngas) * energy

    return pEnergy, cs
