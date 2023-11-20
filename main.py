######################################
# AREPO Initial Conditions Generator #
# mc 18/10/2023 based on work by pcc #
######################################

# Imports
import numpy as np
import os

# Reading in the config and physical properties files
directory = os.path.dirname(os.path.realpath(__file__))

# Getting the path to the config and parameter files
paramFile = str(directory) + "/params.txt"
configFile = str(directory) + "/config.txt"

# Opening the config and parameter files
params = np.loadtxt(paramFile, dtype=float, converters=float)
config = np.loadtxt(configFile, dtype=str)

# Unpacking the parameter file
ngas = int(params[0])
bounds = [params[1], params[2], params[3], params[4], params[5], params[6]]
radius = params[7]
totalMass = params[8]
temperature = params[9]
mu = params[10]
epsilon = params[11]
beta = params[12]
boxDims = [params[13], params[14], params[15]]
tempFactor = params[16]
densityTarget = params[17]

# Defining the code units
uMass = 1.991e33    # grams
uDist = 1e17        # cm
uVelo = 36447.2682  # cm/s
uEner = 1.328e9     # ergs

#######################
# Grid type selection #
#######################

# Uniform particle grid setups
if config[0] == "boxGrid":
    from boxCreation import boxGrid

    # Creating a box grid
    pos, ngas, bounds, volume = boxGrid(ngas, bounds)

    # Running the spherical cut module if sphere selected
elif config[0] == "sphereGrid":
    # Import modules for the box and then spherical grid
    from boxCreation import boxGrid
    from shapeTypes import sphericalCloud

    # Increase ngas as we will lose some particles when cutting out the sphere
    ngas = int(ngas * 6 / np.pi)

    # Creating a box grid
    pos, ngas, bounds, volume = boxGrid(ngas, bounds)

    # Creating spherical grid
    ngas, pos, volume = sphericalCloud(pos, radius, ngas, bounds)

# Randomly placed particle setups
elif config[0] == "boxRan":
    from boxCreation import boxRandom

    # Creating random box grid
    pos, volume = boxRandom(ngas, bounds)

elif config[0] == "sphereRan":
    from boxCreation import sphereRandom

    # Creating a random spherical grid
    pos, volume = sphereRandom(ngas, radius)

# Adjusting positions to be in cm
pos = pos * 3.09e18

###########################
# Mass and energy defines #
###########################

from massAndEnergy import masses
from massAndEnergy import thermalEnergy

# Setting equal particle masses
pMass = masses(ngas, totalMass)

# Converting mass into grams
pMass = pMass * 1.991e33 

# Working out internal energy of each particle along with the sound speed
pEnergy, cs = thermalEnergy(ngas, temperature, mu)

# Converting energy into ergs 
pEnergy = pEnergy * 1e7

################################
# Velocities: Turbulence setup #
################################

# Setup for turbulence from a velocity cube file
if config[1] == "turbFile":
    print("Assigning turbulent velocities")
    from turbulence import turbulenceFromFile

    # Loading in the turbulent velocities from the velocity cube
    velx, vely, velz = turbulenceFromFile(int(config[3]), config[2])

    # Branch for the box scenarios
    if config[0] == "boxGrid" or config[0] == "boxRan":
        from turbulence import boxGridTurbulence

        # Interpolating and assignning velocities
        vels = boxGridTurbulence(velx, vely, velz, pos, pMass, int(config[3]), epsilon)

    # Branch for the spherical scenarios
    elif config[0] == "sphereGrid" or config[0] == "sphereRan":
        from turbulence import sphericalGridTurbulence

        # Interpolating and assigning velocities
        vels = sphericalGridTurbulence(velx, vely, velz, pos, pMass, int(config[3]), epsilon)
    else:
        pass
else:
    # Assgining an empty velocity array if no tubulence setup
    vels = np.zeros((3, ngas), dtype=np.float64)

########################
# Velocities: Rotation #
########################

# Add rotation to the body
if config[4] == "rotation":
    print("Adding solid body rotation")
    from rotation import addRotation

    # Add rotation around z axis of given beta energy ratio
    vels = addRotation(pos, pMass, vels, beta)
else:
    pass

###################################
# Setting particle identification #
###################################

# Assigning each particle an ID from 1 to the max number of particles
pIDs = np.linspace(1, ngas, ngas, dtype=np.int32)

################################
# Low density particle padding #
################################

# Setup for padding the box with low density particles
if config[5] == "True":
    # Branch for the box setups
    print("Padding box with low density particles")
    if config[0] == "boxGrid" or config[0] == "boxRan":
        from lowDensityPadding import padBox

        # Pad the box with low density particles outside the box grid
        pos, vels, pMass, pIDs, pEnergy, pRho, ngasAll = padBox(ngas, pos, vels, pMass, pIDs, pEnergy, boxDims, tempFactor)
    
    # Branch for the spherical setups
    elif config[0] == "sphereGrid" or config[0] == "sphereRan":
        from lowDensityPadding import padSphere

        # Pad the box with low density particles outside the spherical cloud
        pos, vels, pMass, pIDs, pEnergy, pRho, ngasAll = padSphere(ngas, pos, vels, pMass, pIDs, pEnergy, boxDims, tempFactor)
    else:
        ngasAll = ngas
else:
    ngasAll = ngas
    pass

#############################
# Making positions positive #
#############################

# Getting the minimum value of every coordinate
minx = np.min(pos[0])
miny = np.min(pos[1])
minz = np.min(pos[2])

# Shifting everything if its less than zero
if minx < 0:
    pos[0] -= minx
if miny < 0:
    pos[1] -= miny
if minz < 0:
    pos[2] -= minz

############################################
# Conversion of quantities into code units #
############################################

# All variables should be in c.g.s units for conversion
pos = pos / uDist
vels = vels / uVelo
pMass = pMass / uMass
pEnergy = pEnergy / uEner

##############################
# Desired Density Conversion #
##############################

# Converting the number density to code units
densityTarget = densityTarget * 1.4 * 1.66e-24 
densityTarget = densityTarget / (uMass / (uDist**3))
densityTargetPadding = densityTarget * 0.01

# Creating density array
pDensity = np.ones_like(pMass)
pDensity[0:ngas] = densityTarget
pDensity[ngas:-1] = densityTargetPadding

########################
# File output to AREPO #
########################

print("Writing output file")

if config[6] == "hdf5":
    from arepoOut import hdf5out

    # Writing masses to mass 
    if config[7] == "masses":
        # Write the particle data as a hdf5 file
        hdf5out(ngasAll, pos, vels, pIDs, pMass, pEnergy)
    # Writing density to mass
    elif config[7] == "density":
        # Writing particle data as a hdf5 file
        hdf5out(ngasAll, pos, vels, pIDs, pMass, pEnergy, True, pDensity)
else:
    print("Fortran binary version is broken, sorry </3")