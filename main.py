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

paramFile = str(directory) + "\\params.txt"
configFile = str(directory) + "\\config.txt"

params = np.loadtxt(paramFile, dtype=int)
config = np.loadtxt(configFile, dtype=str)

# Unpacking the parameter file
ngas = params[0]
bounds = [params[1], params[2], params[3], params[4], params[5], params[6]]
radius = params[7]
totalMass = params[8]
temperature = params[9]
mu = params[10]
epsilon = params[11]
boxDims = [params[12], params[13], params[14]]
tempFactor = params[15]

#######################
# Grid type selection #
#######################

# Tracking how long each section takes to run
initTime = time()

# Uniform particle grid setups
if config[0] == "boxGrid" or config[0] == "sphereGrid":
    from boxCreation import boxGrid

    # Creating a box grid
    pos, ngas, bounds, volume = boxGrid(ngas, bounds)

    # Running the spherical cut module if sphere selected
    if config[0] == "sphereGrid":
        from shapeTypes import sphericalCloud

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

###########################
# Mass and energy defines #
###########################

from massAndEnergy import masses
from massAndEnergy import thermalEnergy

# Setting equal particle masses
pMass = masses(ngas, totalMass)

# Working out internal energy of each particle along with the sound speed
pEnergy, cs = thermalEnergy(ngas, temperature, mu)

# Time taken for the energy and mass assignment
mAndeTime = time()
print("Particle masses and energies calculated in {:.2f}".format(mAndeTime-gridFinTime))

#################################
# Velocity and Turbulence setup #
#################################

# Setup for turbulence from a velocity cube file
if config[1] == "turbFile":
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
    vels = np.zeros((3, ngas))

# Time taken for the turbulence/velocity setup
turbTime = time()
print("Turbulent velocities assigned in {:.2f}".format(turbTime-mAndeTime))

###################################
# Setting particle identification #
###################################

# ssigning each particle an ID from 0 to the max number of particles
pIDs = np.linspace(0, ngas, ngas, dtype=int)

################################
# Low density particle padding #
################################

# Setup for padding the box with low density particles
if config[4] == "True":
    if config[0] == "boxGrid":
        from lowDensityPadding import padBox

        # Pad the box with low density particles outside the box grid
        pos, vels, pMass, pIDs, pEnergy, pRho = padBox(ngas, pos, vels, pMass, pIDs, pEnergy, boxDims, tempFactor)
    
    elif config[0] == "sphereGrid":
        from lowDensityPadding import padSphere

        # Pad the box with low density particles outside the spherical cloud
        pos, vels, pMass, pIDs, pEnergy, pRho = padSphere(ngas, pos, vels, pMass, pIDs, pEnergy, tempFactor)
else:
    pass

# Time taken to pad the box 
padTime = time()
print("Box padded with low density particles in {:.2f} s".format(padTime-turbTime))

########################
# File output to AREPO #
########################




###########
# TESTING #
###########

import matplotlib.pyplot as plt

''''
x = pos[0,0:ngas]
x2 = pos[0,ngas:]
y = pos[1,0:ngas]
y2 = pos[1,ngas:]
z = pos[2,0:ngas]
z2 = pos[2,ngas:]

fig = plt.figure(figsize=(10,10))
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
x = x[z==np.min(z)]

y = pos[1][:]
y = y[z==np.min(z)]

print(len(x),len(x), np.shape(vels))

fig = plt.figure()
#ax = fig.add_subplot(projection="3d")
ax = fig.add_subplot()

vx = vels[0]
vx = vx[z==np.min(z)]

vy = vels[1]
vy = vy[z==np.min(z)]

ax.quiver(x,y,vx,vy)

#ax.scatter(x, y, z, s=vels[0,:]*500)

#plt.show()

#x = np.linspace(0, 128, 128, endpoint=True)
#y = np.linspace(0, 128, 128, endpoint=True)

#X, Y = np.meshgrid(x, y)

#vx = velx[:,:,0]
#vy = vely[:,:,0]

#ax.quiver(X,Y,vx,vy)
plt.show()
