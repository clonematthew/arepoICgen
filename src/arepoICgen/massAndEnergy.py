# Importing libraries
import numpy as np

# Constants
gasConstant = 8.31

# Defining mass of the particles 
def masses(ngas, config, params):
    # Calculate mass of each particle 
    particleMass = params["mass"] / ngas

    # Set array of particle masses
    pMass = np.ones(ngas, dtype=np.float64) * particleMass

    if config["verbose"]:
        # Printing the mass
        print("Total desired mass: %s" % params["mass"])
        print("Initial particle mass: {:.5f}".format(particleMass))

    return pMass

# Defining energy of the particles
def thermalEnergy(ngas, config, params):
    # Calculating internal energy
    particleMass = params["mass"] * 1.991e33 / ngas
    
    energy = (3./2.) * particleMass * params["temp"] * gasConstant / params["mu"]

    # Calculating a sound speed
    cs = np.sqrt(gasConstant * params["temp"])

    # Allocating this energy to each particle
    pEnergy = np.ones(ngas, dtype=np.float64) * energy

    if config["verbose"]:
        # Printing values
        print("Initial particle energy: {:.2e} ergs".format(energy))
        print("Initial sound speed: {:.2f} m/s".format(cs))

    return pEnergy
