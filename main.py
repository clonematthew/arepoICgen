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

#######################
# Grid type selection #
#######################

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
        ngas, pos, volume = sphericalCloud(pos, radius)

# Randomly placed particle setups
elif config[0] == "boxRan":
    from boxCreation import boxRandom

    # Creating random box grid
    pos, volume = boxRandom(ngas, bounds)

elif config[0] == "sphereRan":
    from boxCreation import sphereRandom

    # Creating a random spherical grid
    pos, volume = sphereRandom(ngas, radius)

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

    if config[0] == "boxGrid" or config[0] == "sphereGrid":
        from turbulence import boxGridTurbulence

        # Interpolating and assignning velocities
        vels = boxGridTurbulence(velx, vely, velz, pos, pMass, int(config[3]), epsilon)
    
    else:
        pass

else:
    # Assgining an empty velocity array if no tubulence setup
    vels = np.zeros((3, ngas))

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
    if config[0] == "sphereGrid":
        from lowDensityPadding import padBoxSphere

        # Pad the box with low density particles outside the spherical cloud
        pos, vels, pMass, pIDs, pEnergy, pRho = padBoxSphere(ngas, pos, vels, pMass, pIDs, pEnergy)
    
    else:
        pass
else:
    pass

###########
# TESTING #
###########

import matplotlib.pyplot as plt

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