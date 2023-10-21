# Needed libraries
import numpy as np
from random import random

# Code to create particles in a box grid
def boxGrid(ngas, bounds):
    # Unpacking bounds
    xmin = bounds[0]
    xmax = bounds[1]
    ymin = bounds[2]
    ymax = bounds[3]
    zmin = bounds[4]
    zmax = bounds[5]

    # Calculating the box volume
    volume = (xmax-xmin) * (ymax-ymin) * (zmax-zmin)

    # Determining average particle spacing
    spacing = (volume / ngas)**(1./3.)

    # Finding the number of grid points for each dimension
    nx = int((xmax-xmin)/spacing)
    ny = int((ymax-ymin)/spacing)
    nz = int((zmax-zmin)/spacing)

    # Resetting the number of gas particles to the rounded version 
    ngas = nx * ny * nz

    # Creating arrays for the particles
    pos = np.zeros((3, ngas))

    # Looping through every particle and assigning its position
    npart = 0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                # Assigning positions 
                pos[0][npart] = xmin + (0.5 * spacing) + (i * spacing)
                pos[1][npart] = ymin + (0.5 * spacing) + (j * spacing)
                pos[2][npart] = zmin + (0.5 * spacing) + (k * spacing)

                # Adding to counter
                npart += 1

    # Setting the max dimensions to the maximum particle positions 
    xmin = min(pos[0][:])
    xmax = max(pos[0][:])
    ymin = min(pos[1][:])
    ymax = max(pos[1][:])
    zmin = min(pos[2][:])
    zmax = max(pos[2][:]) 

    # Calculating an updated volume using these
    volume = (xmax-xmin) * (ymax-ymin) * (zmax-zmin)

    return pos, ngas, [xmin, xmax, ymin, ymax, zmin, zmax], volume

# Code to allocate particle positions randomly
def boxRandom(ngas, bounds):
    # Unpacking bounds
    xmin = bounds[0]
    xmax = bounds[1]
    ymin = bounds[2]
    ymax = bounds[3]
    zmin = bounds[4]
    zmax = bounds[5]
    
    # Creating the particle array
    pos = np.zeros((3, ngas))

    # Calculating volume
    volume = (xmax-xmin) * (ymax-ymin) * (zmax-zmin)

    # Looping through the list of particles
    for i in range(ngas):
        pos[0][i] = xmin + (xmax - xmin) * random()
        pos[1][i] = ymin + (ymax - ymin) * random()
        pos[2][i] = zmin + (zmax - zmin) * random()

    return pos, volume

def sphereRandom(ngas, radius):
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

        if r <= radius:
            pos[0][i] = x
            pos[1][i] = y
            pos[2][i] = z

            i += 1

    return pos, volume


