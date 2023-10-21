# Importing libraries
import numpy as np
from scipy import interpolate

# Function to load in the velocities from file
def turbulenceFromFile(gridSize, filename):
    # Calculating number of cells
    nCells = gridSize - 1 

    # Allocating arrays for the velocities
    velx = np.zeros((3, gridSize))
    vely = np.zeros((3, gridSize))
    velz = np.zeros((3, gridSize))

    # Opening the file containing the velocity cube
    with open(filename, "r") as f:
        # Set vel x
        # Set vel y
        # Set vel z
        test = 1    

    return velx, vely, velz

# Function to interpolate velocities for a box grid
def boxGridTurbulence(velx, vely, velz, pos, pMass, gridSize, epsilon):
    # Setting up centre of mass positions
    xcom = 0.
    ycom = 0.
    zcom = 0.
    mtot = 0.

    # Finding the centre of mass
    ngas = len(pos[0][:])
    for i in range(ngas):
        # Calculating centre of masses
        mtot += pMass[i]
        xcom += pos[0][i] * pMass[i]
        ycom += pos[1][i] * pMass[i]
        zcom += pos[2][i] * pMass[i]

    # Scaling centre of mass
    xcom = xcom / mtot
    ycom = ycom / mtot
    zcom = zcom / mtot

    # Finding radial positions of each particle w/r/t to the centre of mass
    for i in range(ngas):
        # Calculating r
        r = (pos[0][i] - xcom)**2 + (pos[1][i] - ycom)**2 + (pos[2][i] - zcom)**2

    # Setting radnorm to the further away particle from CoM
    radnorm = max(r)

    deli = radnorm / gridSize  

    # Creating the boxGrid for the velocities 
    x = np.linspace(min(pos[0][:]), max(pos[0][:]), gridSize)
    y = np.linspace(min(pos[1][:]), max(pos[1][:]), gridSize)
    z = np.linspace(min(pos[2][:]), max(pos[2][:]), gridSize)

    xg, yg, zg = np.meshgrid(x, y, z)

    # Setting up interpolation objects for each dimension
    xInterp = interpolate.RegularGridInterpolator((xg, yg, zg), velx)
    yInterp = interpolate.RegularGridInterpolator((xg, yg, zg), vely)
    zInterp = interpolate.RegularGridInterpolator((xg, yg, zg), velz)

    # Creating velocity arrays for the particles
    vels = np.zeros((3, ngas))
    velxcom = 0.
    velycom = 0.
    velzcom = 0.

    # Assigning a velocity for each particle via interpolation of the grid
    for i in range(ngas):
        # Getting current position
        currentPos = np.array([pos[0][i], pos[1][i], pos[2][i]])

        # Interpolating the velocities in each dimension
        xvel = xInterp(currentPos)
        yvel = yInterp(currentPos)
        zvel = zInterp(currentPos)

        # Applying velocities to the velocity array
        vels[0][i] = xvel
        vels[1][i] = yvel
        vels[2][i] = zvel

        # Adding to centre of mass velocity values
        velxcom += xvel * pMass[i]
        velycom += yvel * pMass[i]
        velzcom += zvel * pMass[i]

    # Finding the true centre of mass velocity 
    velxcom = velxcom / mtot
    velycom = velycom / mtot
    velzcom = velzcom / mtot

    # Subtracting the velocity from each particle
    for i in range(ngas):
        vels[0][i] -= velxcom
        vels[1][i] -= velycom
        vels[2][i] -= velzcom

    # Finding the total kinetic energy of the cloud
    eKinetic = 0.
    rootMeanSquared = 0.
    vmax = 0.

    for i in range(ngas):
        # Find current particles velocity and add to counters
        velAbs = vels[0][i]**2 + vels[1][i]**2 + vels[2][i]**2
        rootMeanSquared += velAbs
        eKinetic += velAbs*pMass[i]*0.5

        # Update the maxium velocity
        vmax = max([vmax, vels[0][i], vels[1][i], vels[2][i]])

    # Calculate the gravitational potential of the cloud
    ePotential = (3/5) * mtot * mtot / radnorm

    # Calculate the root mean squared speed of the cloud
    rootMeanSquared = np.sqrt(rootMeanSquared / ngas)

    ## FOR SCALING BY VIRIAL EQUILIBRIUM
    scalingFactor = (ePotential/eKinetic) / epsilon

    # Scaling velocities
    vels = vels * scalingFactor

    return vels

    








    