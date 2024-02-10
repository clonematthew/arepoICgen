# Importing libraries
import numpy as np
from scipy.io import FortranFile
from numba import jit

# Function to load in the velocities from file
def turbulenceFromFile(gridSize, filename):
    # Loading in the fortran file using sciPy
    f = FortranFile(filename, "r")

    # Extracting each block and assignning to arrays
    velx = f.read_reals().reshape((gridSize,gridSize,gridSize))
    vely = f.read_reals().reshape((gridSize,gridSize,gridSize))
    velz = f.read_reals().reshape((gridSize,gridSize,gridSize))

    return velx, vely, velz

# Function to carry out the interpolation, compiled with jit
@jit(nopython=True)
def turbLoopBox(ngas, pos, vels, radnorm, gridSize, velx, vely, velz, deli):
    for i in range(ngas):
        # Finding positions in the velocity array our particle positions refer to
        iposx = np.int64((pos[0,i]/radnorm)*gridSize+0.5)
        iposx = np.min(np.array([np.max(np.array([iposx, 0.], dtype=np.int64)), gridSize-2.], dtype=np.int64))

        iposy = np.int64((pos[1,i]/radnorm)*gridSize+0.5)
        iposy = np.min(np.array([np.max(np.array([iposy, 0.], dtype=np.int64)), gridSize-2.], dtype=np.int64))

        iposz = np.int64((pos[2,i]/radnorm)*gridSize+0.5)
        iposz = np.min(np.array([np.max(np.array([iposz, 0.], dtype=np.int64)), gridSize-2.], dtype=np.int64))
    
        delx = pos[0,i] - (iposx - 0.5) / np.real(gridSize) * radnorm
        dely = pos[1,i] - (iposy - 0.5) / np.real(gridSize) * radnorm
        delz = pos[2,i] - (iposz - 0.5) / np.real(gridSize) * radnorm

        # Interpolating the x velocity
        velx1 = velx[iposx,iposy,iposz] + delx/deli * ((velx[iposx+1,iposy,iposz]) - velx[iposx,iposy,iposz])
        velx2 = velx[iposx,iposy+1,iposz] + delx/deli * (velx[iposx+1,iposy+1,iposz] - velx[iposx,iposy+1,iposz])
        vely1 = velx1 + dely/deli*(velx2-velx1)

        velx1 = velx[iposx,iposy,iposz+1] + delx/deli * (velx[iposx+1,iposy,iposz+1]-velx[iposx,iposy,iposz+1])
        velx2 = velx[iposx,iposy+1,iposz+1] + delx/deli * (velx[iposx+1,iposy+1,iposz+1]-velx[iposx,iposy+1,iposz+1])
        vely2 = velx1 + dely/deli*(velx2-velx1)

        vels[0][i] = vely1 + delz/deli*(vely2-vely1)

        # Interpolating the y velocity
        velx1 = vely[iposx,iposy,iposz] + delx/deli * ((vely[iposx+1,iposy, iposz]) - vely[iposx,iposy,iposz])
        velx2 = vely[iposx,iposy+1,iposz] + delx/deli * (vely[iposx+1,iposy+1,iposz] - vely[iposx,iposy+1,iposz])
        vely1 = velx1 + dely/deli*(velx2-velx1)

        velx1 = vely[iposx,iposy,iposz+1] + delx/deli * (vely[iposx+1,iposy,iposz+1]-vely[iposx,iposy,iposz+1])
        velx2 = vely[iposx,iposy+1,iposz+1] + delx/deli * (vely[iposx+1,iposy+1,iposz+1]-vely[iposx,iposy+1,iposz+1])
        vely2 = velx1 + dely/deli*(velx2-velx1)

        vels[1][i] = vely1 + delz/deli*(vely2-vely1)

        # Interpolating the z velocity
        velx1 = velz[iposx,iposy,iposz] + delx/deli * ((velz[iposx+1,iposy, iposz]) - velz[iposx,iposy,iposz])
        velx2 = velz[iposx,iposy+1,iposz] + delx/deli * (velz[iposx+1,iposy+1,iposz] - velz[iposx,iposy+1,iposz])
        vely1 = velx1 + dely/deli*(velx2-velx1)

        velx1 = velz[iposx,iposy,iposz+1] + delx/deli * (velz[iposx+1,iposy,iposz+1]-velz[iposx,iposy,iposz+1])
        velx2 = velz[iposx,iposy+1,iposz+1] + delx/deli * (velz[iposx+1,iposy+1,iposz+1]-velz[iposx,iposy+1,iposz+1])
        vely2 = velx1 + dely/deli*(velx2-velx1)

        vels[2][i] = vely1 + delz/deli*(vely2-vely1)

    return vels

# Function to interpolate velocites for the spherical grid
@jit(nopython=True)
def turbLoopSphere(ngas, pos, vels, radnorm, gridSize, velx, vely, velz, deli):
    for i in range(ngas):
        # Finding positions in the velocity array our particle positions refer to
        iposx = np.int64(pos[0,i]/radnorm * (gridSize/2.) + (gridSize/2.) + 0.5)
        iposx = np.min(np.array([np.max(np.array([iposx, 0.], dtype=np.int64)), gridSize-2.], dtype=np.int64))

        iposy = np.int64(pos[1,i]/radnorm * (gridSize/2.) + (gridSize/2.) + 0.5)
        iposy = np.min(np.array([np.max(np.array([iposy, 0.], dtype=np.int64)), gridSize-2.], dtype=np.int64))

        iposz = np.int64(pos[2,i]/radnorm * (gridSize/2.) + (gridSize/2.) + 0.5)
        iposz = np.min(np.array([np.max(np.array([iposz, 0.], dtype=np.int64)), gridSize-2.], dtype=np.int64))

        delx = pos[0,i] - (iposx - (gridSize/2.) - 0.5) / np.real(gridSize/2.) * radnorm
        dely = pos[1,i] - (iposy - (gridSize/2.) - 0.5) / np.real(gridSize/2.) * radnorm
        delz = pos[2,i] - (iposz - (gridSize/2.) - 0.5) / np.real(gridSize/2.) * radnorm

        # Interpolating the x velocity
        velx1 = velx[iposx,iposy,iposz] + delx/deli * ((velx[iposx+1,iposy,iposz]) - velx[iposx,iposy,iposz])
        velx2 = velx[iposx,iposy+1,iposz] + delx/deli * (velx[iposx+1,iposy+1,iposz] - velx[iposx,iposy+1,iposz])
        vely1 = velx1 + dely/deli*(velx2-velx1)

        velx1 = velx[iposx,iposy,iposz+1] + delx/deli * (velx[iposx+1,iposy,iposz+1]-velx[iposx,iposy,iposz+1])
        velx2 = velx[iposx,iposy+1,iposz+1] + delx/deli * (velx[iposx+1,iposy+1,iposz+1]-velx[iposx,iposy+1,iposz+1])
        vely2 = velx1 + dely/deli*(velx2-velx1)

        vels[0][i] = vely1 + delz/deli*(vely2-vely1)

        # Interpolating the y velocity
        velx1 = vely[iposx,iposy,iposz] + delx/deli * ((vely[iposx+1,iposy, iposz]) - vely[iposx,iposy,iposz])
        velx2 = vely[iposx,iposy+1,iposz] + delx/deli * (vely[iposx+1,iposy+1,iposz] - vely[iposx,iposy+1,iposz])
        vely1 = velx1 + dely/deli*(velx2-velx1)

        velx1 = vely[iposx,iposy,iposz+1] + delx/deli * (vely[iposx+1,iposy,iposz+1]-vely[iposx,iposy,iposz+1])
        velx2 = vely[iposx,iposy+1,iposz+1] + delx/deli * (vely[iposx+1,iposy+1,iposz+1]-vely[iposx,iposy+1,iposz+1])
        vely2 = velx1 + dely/deli*(velx2-velx1)

        vels[1][i] = vely1 + delz/deli*(vely2-vely1)

        # Interpolating the z velocity
        velx1 = velz[iposx,iposy,iposz] + delx/deli * ((velz[iposx+1,iposy, iposz]) - velz[iposx,iposy,iposz])
        velx2 = velz[iposx,iposy+1,iposz] + delx/deli * (velz[iposx+1,iposy+1,iposz] - velz[iposx,iposy+1,iposz])
        vely1 = velx1 + dely/deli*(velx2-velx1)

        velx1 = velz[iposx,iposy,iposz+1] + delx/deli * (velz[iposx+1,iposy,iposz+1]-velz[iposx,iposy,iposz+1])
        velx2 = velz[iposx,iposy+1,iposz+1] + delx/deli * (velz[iposx+1,iposy+1,iposz+1]-velz[iposx,iposy+1,iposz+1])
        vely2 = velx1 + dely/deli*(velx2-velx1)

        vels[2][i] = vely1 + delz/deli*(vely2-vely1)

    return vels  

# Function to centre and scale velocities
def scaleVelocities(ngas, vels, pMass, radnorm, epsilon):
    # Finding total mass
    mtot = np.sum(pMass)

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
    absVel = (vels[0]**2 + vels[1]**2 + vels[2]**2)
    eKinetic = np.sum(absVel * pMass * 0.5)

    # Calculate the gravitational potential of the cloud
    ePotential = 6.67e-8 * (3/5) * mtot * mtot / radnorm

    ## FOR SCALING BY VIRIAL EQUILIBRIUM
    scalingFactor = np.sqrt((ePotential/eKinetic) / epsilon)

    # Scaling velocities
    vels = vels * scalingFactor

    return vels

# Function to interpolate velocities for a box grid
def boxGridTurbulence(velx, vely, velz, pos, pMass, gridSize, epsilon):
    # Moving the box to have no negative values
    xmin = np.min(pos[0])
    if xmin < 0:
        pos[0] += abs(xmin)

    ymin = np.min(pos[1])
    if ymin < 0:
        pos[1] += abs(ymin)
    
    zmin = np.min(pos[2])
    if zmin < 0:
        pos[2] += abs(zmin)

    # Finding the centre of mass
    ngas = len(pMass)

    # Finding scaling factors
    radnorm = np.max([np.max(pos[0]), np.max(pos[1]), np.max(pos[2])])
    deli = radnorm / gridSize  

    # Creating velocity arrays for the particles
    vels = np.zeros((3, ngas), dtype=np.float64)

    # Interpolating the velocities ## BOX GRID METHOD
    vels = turbLoopBox(ngas, pos, vels, radnorm, gridSize, velx, vely, velz, deli)

    # Scaling the velocities by the required amount
    vels = scaleVelocities(ngas, vels, pMass, radnorm, epsilon)

    return vels

# Function to interpolate the velocities for a spherical grid
def sphericalGridTurbulence(velx, vely, velz, pos, pMass, gridSize, epsilon):

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
    radnorm = np.sqrt(np.max(r))
    deli = radnorm/(gridSize/2)

    # Creating velocity arrays for the particles
    vels = np.zeros((3, ngas), dtype=np.double)

    # Interpolating velocities with the spherical method
    vels = turbLoopSphere(ngas, pos, vels, radnorm, gridSize, velx, vely, velz, deli)
    
    # Scaling the velocities by the required amount
    vels = scaleVelocities(ngas, vels, pMass, radnorm, epsilon)

    return vels



    








    