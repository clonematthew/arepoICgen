# Importing libraries
import numpy as np
from scipy import interpolate
from scipy.io import FortranFile

# Function to load in the velocities from file
def turbulenceFromFile(gridSize, filename):
    # Calculating number of cells
    nCells = gridSize - 1 

    # Loading in the fortran file using sciPy
    f = FortranFile(filename, "r")

    # Extracting each block and assignning to arrays
    velx = f.read_reals().reshape((gridSize,gridSize,gridSize))
    vely = f.read_reals().reshape((gridSize,gridSize,gridSize))
    velz = f.read_reals().reshape((gridSize,gridSize,gridSize))

    return velx, vely, velz

# Function to interpolate velocities for a box grid
def boxGridTurbulence(velx, vely, velz, pos, pMass, gridSize, epsilon):
    # Setting up centre of mass positions
    xcom = 0.
    ycom = 0.
    zcom = 0.
    mtot = 0.

    # Finding the centre of mass
    ngas = len(pMass)

    mtot = np.sum(pMass)
    xcom = np.sum(pos[0]*pMass)
    ycom = np.sum(pos[1]*pMass)
    zcom = np.sum(pos[2]*pMass)

    # Scaling centre of mass
    xcom = xcom / mtot
    ycom = ycom / mtot
    zcom = zcom / mtot

    # Finding radial positions of each particle w/r/t to the centre of mass
    r = (pos[0] - xcom)**2 + (pos[1] - ycom)**2 + (pos[2] - zcom)**2

    # Setting radnorm to the further away particle from CoM
    radnorm = np.max(r)

    deli = radnorm / gridSize  


    # Creating velocity arrays for the particles
    vels = np.zeros((3, ngas))

    '''
    # Creating the boxGrid for the velocities 
    x = np.linspace(min(pos[0][:]), max(pos[0][:]), gridSize)
    y = np.linspace(min(pos[1][:]), max(pos[1][:]), gridSize)
    z = np.linspace(min(pos[2][:]), max(pos[2][:]), gridSize)

    # Setting up interpolation objects for each dimension
    xInterp = interpolate.RegularGridInterpolator(points=[x,y,z], values=velx)
    yInterp = interpolate.RegularGridInterpolator(points=[x,y,z], values=vely)
    zInterp = interpolate.RegularGridInterpolator(points=[x,y,z], values=velz)


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

    '''

    # Trying manual interpolation
    for i in range(ngas):
        iposx = int(pos[0][i]/radnorm*(gridSize/2)+(gridSize/2)+0.5)
        iposx = np.min([np.max([iposx,1]), gridSize-1])

        iposy = int(pos[1][i]/radnorm*(gridSize/2)+(gridSize/2)+0.5)
        iposy = np.min([np.max([iposy,1]), gridSize-1])

        iposz = int(pos[2][i]/radnorm*(gridSize/2)+(gridSize/2)+0.5)
        iposz = np.min([np.max([iposz,1]), gridSize-1])

        delx = pos[0][i] - (iposx - (gridSize/2) - 0.5) / np.real(gridSize/2)*radnorm
        dely = pos[1][i] - (iposy - (gridSize/2) - 0.5) / np.real(gridSize/2)*radnorm
        delz = pos[2][i] - (iposz - (gridSize/2) - 0.5) / np.real(gridSize/2)*radnorm
        
        # x vel
        velx1 = velx[iposx,iposy,iposz] + delx/deli * ((velx[iposx+1,iposy,iposz]) - velx[iposx,iposy,iposz])
        velx2 = velx[iposx,iposy+1,iposz] + delx/deli * (velx[iposx+1,iposy+1,iposz] - velx[iposx,iposy+1,iposz])
        vely1 = velx1 + dely/deli*(velx2-velx1)

        velx1 = velx[iposx,iposy,iposz+1] + delx/deli * (velx[iposx+1,iposy,iposz+1]-velx[iposx,iposy,iposz+1])
        velx2 = velx[iposx,iposy+1,iposz+1] + delx/deli * (velx[iposx+1,iposy+1,iposz+1]-velx[iposx,iposy+1,iposz+1])
        vely2 = velx1 + dely/deli*(velx2-velx1)

        vels[0][i] = vely1 + delz/deli*(vely2-vely1)

        # y vel
        velx1 = vely[iposx,iposy,iposz] + delx/deli * ((vely[iposx+1,iposy, iposz]) - vely[iposx,iposy,iposz])
        velx2 = vely[iposx,iposy+1,iposz] + delx/deli * (vely[iposx+1,iposy+1,iposz] - vely[iposx,iposy+1,iposz])
        vely1 = velx1 + dely/deli*(velx2-velx1)

        velx1 = vely[iposx,iposy,iposz+1] + delx/deli * (vely[iposx+1,iposy,iposz+1]-vely[iposx,iposy,iposz+1])
        velx2 = vely[iposx,iposy+1,iposz+1] + delx/deli * (vely[iposx+1,iposy+1,iposz+1]-vely[iposx,iposy+1,iposz+1])
        vely2 = velx1 + dely/deli*(velx2-velx1)

        vels[1][i] = vely1 + delz/deli*(vely2-vely1)

        # z vel
        velx1 = velz[iposx,iposy,iposz] + delx/deli * ((velz[iposx+1,iposy, iposz]) - velz[iposx,iposy,iposz])
        velx2 = velz[iposx,iposy+1,iposz] + delx/deli * (velz[iposx+1,iposy+1,iposz] - velz[iposx,iposy+1,iposz])
        vely1 = velx1 + dely/deli*(velx2-velx1)

        velx1 = velz[iposx,iposy,iposz+1] + delx/deli * (velz[iposx+1,iposy,iposz+1]-velz[iposx,iposy,iposz+1])
        velx2 = velz[iposx,iposy+1,iposz+1] + delx/deli * (velz[iposx+1,iposy+1,iposz+1]-velz[iposx,iposy+1,iposz+1])
        vely2 = velx1 + dely/deli*(velx2-velx1)

        vels[2][i] = vely1 + delz/deli*(vely2-vely1)

    # Finding centre of mass velocity
    velxcom = np.sum(vels[0]*pMass)
    velycom = np.sum(vels[1]*pMass)
    velzcom = np.sum(vels[2]*pMass)

    velxcom = velxcom / mtot
    velycom = velycom / mtot
    velzcom = velzcom / mtot

    # Subtracting the velocity from each particle
    vels[0] = vels[0] - velxcom
    vels[1] = vels[1] - velycom
    vels[2] = vels[2] - velzcom

    # Finding the total kinetic energy of the cloud
    absVel = vels[0]**2 + vels[1]**2 + vels[2]**2
    rootMeanSquared = np.sum(absVel)
    eKinetic = np.sum(absVel*pMass*0.5)

    # Calculate the gravitational potential of the cloud
    ePotential = (3/5) * mtot * mtot / radnorm

    # Calculate the root mean squared speed of the cloud
    rootMeanSquared = np.sqrt(rootMeanSquared / ngas)

    ## FOR SCALING BY VIRIAL EQUILIBRIUM
    scalingFactor = (ePotential/eKinetic) / epsilon

    # Scaling velocities
    vels = vels * scalingFactor

    return vels

    








    