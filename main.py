######################################
# AREPO Initial Conditions Generator #
# mc 18/10/2023 based on work by pcc #
######################################

# Imports
import numpy as np
import os
from time import time

# Reading in the config and physical properties files
directory = os.path.dirname(os.path.realpath(__file__))

#paramFile = str(directory) + "\\" + input("Physical Parameter file: ") + ".txt"
#configFile = str(directory) + "\\" + input("Configuration File: ") + ".txt"

paramFile = str(directory) + "/params.txt"
configFile = str(directory) + "/config.txt"

params = np.loadtxt(paramFile, dtype=float, converters=float)
config = np.loadtxt(configFile, dtype=str)

# Unpacking the parameter file
ngas = params[0]
bounds = [params[1], params[2], params[3], params[4], params[5], params[6]]
radius = params[7]
totalMass = params[8]
temperature = params[9]
mu = params[10]
epsilon = params[11]
beta = params[12]
boxDims = [params[13], params[14], params[15]]
tempFactor = params[16]

# Defining the code units
uMass = 1.991e33    # grams
uDist = 1e17        # cm
uVelo = 36447.2682  # cm/s
uEner = 2.64485e42  # ergs

#######################
# Grid type selection #
#######################

# Tracking how long each section takes to run
print("Creating particle grid points")
initTime = time()

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

# Time taken for the grid setup
gridFinTime = time()
print("Particle grid created in {:.2f} s".format(gridFinTime-initTime))

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

# Time taken for the energy and mass assignment
mAndeTime = time()
print("Particle masses and energies calculated in {:.2f}".format(mAndeTime-gridFinTime))

################################
# Velocities: Turbulence setup #
################################

# Setup for turbulence from a velocity cube file
if config[1] == "turbFile":
    print("Assigning turbulent velocities")
    from turbulence import turbulenceFromFile

    # Loading in the turbulent velocities from the velocity cube
    velx, vely, velz = turbulenceFromFile(int(config[3]), config[2])

    if config[0] == "boxGrid":
        from turbulence import boxGridTurbulence

        # Interpolating and assignning velocities
        vels = boxGridTurbulence(velx, vely, velz, pos, pMass, int(config[3]), epsilon)

    elif config[0] == "sphereGrid":
        from turbulence import sphericalGridTurbulence

        # Interpolating and assigning velocities
        vels = sphericalGridTurbulence(velx, vely, velz, pos, pMass, int(config[3]), epsilon)
    
    else:
        pass

else:
    # Assgining an empty velocity array if no tubulence setup
    vels = np.zeros((3, ngas), dtype=np.float64)

# Time taken for the turbulence/velocity setup
turbTime = time()
print("Turbulent velocities assigned in {:.2f}".format(turbTime-mAndeTime))

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

# Assigning each particle an ID from 0 to the max number of particles
pIDs = np.linspace(1, ngas+1, ngas, dtype=np.int32)

################################
# Low density particle padding #
################################

# Setup for padding the box with low density particles
if config[5] == "True":
    print("Padding box with low density particles")
    if config[0] == "boxGrid":
        from lowDensityPadding import padBox

        # Pad the box with low density particles outside the box grid
        pos, vels, pMass, pIDs, pEnergy, pRho, ngasAll = padBox(ngas, pos, vels, pMass, pIDs, pEnergy, boxDims, tempFactor)
    
    elif config[0] == "sphereGrid":
        from lowDensityPadding import padSphere

        # Pad the box with low density particles outside the spherical cloud
        pos, vels, pMass, pIDs, pEnergy, pRho, ngasAll = padSphere(ngas, pos, vels, pMass, pIDs, pEnergy, boxDims, tempFactor)
else:
    pass

# Time taken to pad the box 
padTime = time()
print("Box padded with low density particles in {:.2f} s".format(padTime-turbTime))

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

########################
# File output to AREPO #
########################

print("Writing output file")

if config[6] == "arepo":
    from arepoOut import arepoOut

    # Write out the data to a type 2 AREPO file 
    arepoOut(ngasAll, pos, vels, pIDs, pMass, pEnergy)
else:
    pass

###########
# TESTING #
###########

import matplotlib.pyplot as plt


x = pos[0,0:ngas]
x2 = pos[0,ngas:]
y = pos[1,0:ngas]
y2 = pos[1,ngas:]
z = pos[2,0:ngas]
z2 = pos[2,ngas:]

fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(projection="3d")
ax.scatter(x, y, z, s=1, color="blue")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.scatter(x2, y2, z2, s=1, color="red", )
#plt.scatter(x[z==z[500]], y[z==z[500]], s=1, color="blue")
#plt.scatter(x2[z2==z[500]], y2[z2==z[500]], s=1, color="red")
plt.show()



'''


z = pos[2][:]

x = pos[0][:]
x = x[z==z[20000]]

y = pos[1][:]
y = y[z==z[20000]]

print(len(x),len(x), np.shape(vels))

fig = plt.figure()

#ax = fig.add_subplot(projection="3d")
ax = fig.add_subplot()
ax.set_aspect("equal")

vx = vels[0]
vx = vx[z==z[20000]]

vy = vels[1]
vy = vy[z==z[20000]]

ax.quiver(x,y,vx,vy)
'''

#ax.scatter(x, y, z, s=vels[0,:]*500)

#plt.show()

#x = np.linspace(0, 128, 128, endpoint=True)
#y = np.linspace(0, 128, 128, endpoint=True)

#X, Y = np.meshgrid(x, y)

#vx = velx[:,:,0]
#vy = vely[:,:,0]

#ax.quiver(X,Y,vx,vy)
#plt.show()

