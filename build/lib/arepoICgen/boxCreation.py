# Needed libraries
import numpy as np
from random import random
from numba import jit

# Compiled loop to speed things up, this is slow in normal python for large ngas
@jit(nopython=True)
def tripleLoop(nx, ny, nz, pos, spacing):
    npart = 0
    for i in range(0, nx):
       for j in range(0, ny):
            for k in range(0, nz):
                # Assigning positions 
                pos[0][npart] = (0.5 * spacing) + (i * spacing)
                pos[1][npart] = (0.5 * spacing) + (j * spacing)
                pos[2][npart] = (0.5 * spacing) + (k * spacing)

                # Adding to counter
                npart += 1
    return pos

# Code to create particles in a box grid
def boxGrid(ngas, lengths, verbose=False):
    # Unpacking bounds
    xLength = lengths[0]
    yLength = lengths[1]
    zLength = lengths[2]

    # Calculating the box volume
    volume = xLength * yLength * zLength

    # Determining average particle spacing
    spacing = (volume / ngas)**(1./3.)

    # Finding the number of grid points for each dimension
    nx = np.int64(xLength/spacing)
    ny = np.int64(yLength/spacing)
    nz = np.int64(zLength/spacing)

    # Resetting the number of gas particles to the rounded version 
    ngas = nx * ny * nz

    # Creating arrays for the particles
    pos = np.zeros((3, ngas), dtype=np.float64)

    # Looping through every particle and assigning its position
    pos = tripleLoop(nx, ny, nz, pos, spacing)

    # Setting the max dimensions to the maximum particle positions 
    xmin = np.min(pos[0])
    xmax = np.max(pos[0])
    ymin = np.min(pos[1])
    ymax = np.max(pos[1])
    zmin = np.min(pos[2])
    zmax = np.max(pos[2]) 

    # Calculating an updated volume using these
    volume = (xmax-xmin) * (ymax-ymin) * (zmax-zmin)

    if verbose:
        print("Number of points in each dimension: %s, %s, %s" % (nx, ny, nz))
        print("New number of particles: %s" % ngas)
        print("Spacing between points: {:.2f}".format(spacing))

        print("New X Limits: {:.2f} - {:.2f}".format(xmin, xmax))
        print("New Y Limits: {:.2f} - {:.2f}".format(ymin, ymax))
        print("New Z Limits: {:.2f} - {:.2f}".format(zmin, zmax))

        print("Box Grid Volume: {:.2e}".format(volume))

    return pos, ngas, volume

# Code to allocate particle positions randomly
def boxRandom(ngas, lengths, verbose=False):
    # Unpacking bounds
    xLength = lengths[0]
    yLength = lengths[1]
    zLength = lengths[2]
    
    # Creating the particle array
    pos = np.zeros((3, int(ngas)), dtype=np.float64)

    # Calculating volume
    volume = xLength * yLength * zLength

    # Looping through the list of particles
    for i in range(ngas):
        pos[0][i] = xLength * random()
        pos[1][i] = yLength * random()
        pos[2][i] = zLength * random()

    return pos, volume

def sphereRandom(ngas, radius, verbose=False):
    # Creating particle array
    pos = np.zeros((3, ngas))

    # Calculating volume
    volume = (radius**3) *  (4. * np.pi) / 3.

    # Allocating positions
    i = 0
    while i < ngas:
        x = -radius + 2. * radius * random()    
        y = -radius + 2. * radius * random()
        z = -radius + 2. * radius * random()
        r = np.sqrt(x**2 + y**2 + z**2)

        if x == 0 or y == 0 or z == 0:
            pass
        else:
            if r <= radius:
                pos[0][i] = x
                pos[1][i] = y
                pos[2][i] = z

                i += 1
            else:
                pass

    # Setting the max dimensions to the maximum particle positions 
    xmin = np.min(pos[0])
    xmax = np.max(pos[0])
    ymin = np.min(pos[1])
    ymax = np.max(pos[1])
    zmin = np.min(pos[2])
    zmax = np.max(pos[2]) 

    if verbose:
        # Printing values
        print("Spherical Volume: {:.2e}".format(volume))

        # Printing the new limits
        print("New X Limits: {:.2f} - {:.2f}".format(xmin, xmax))
        print("New Y Limits: {:.2f} - {:.2f}".format(ymin, ymax))
        print("New Z Limits: {:.2f} - {:.2f}".format(zmin, zmax))

    return pos, volume


