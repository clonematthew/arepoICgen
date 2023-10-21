# Needed libraries
import numpy as np

# Function for a spherical cloud
def sphericalCloud(pos, cloudSize):
    # Finding num of particles
    ngas = len(pos[0][:])

    # Getting the min and max positions
    xmin = min(pos[0][:])
    xmax = max(pos[0][:])
    ymin = min(pos[1][:])
    ymax = max(pos[1][:])
    zmin = min(pos[2][:])
    zmax = max(pos[2][:])

    # Getting box sizes
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin

    # Finding the radius of the central sphere
    radius = min([dx/2, dy/2, dz/2])

    # Fiding the central point
    xcom = xmin + dx/2
    ycom = ymin + dy/2
    zcom = zmin + dz/2

    # Only keeping the particles within the spherical region
    newPos = np.zeros((3, ngas))
    ngasNew = 0

    for i in range(ngas):
        # Distance of each particle from centre
        r = (pos[0][i] - xcom)**2 + (pos[1][i] - ycom)**2 + (pos[2][i] - zcom)**2

        # Check if within sphere
        if r <= radius**2:
            # Update the number of particles inside the sphere
            ngasNew += 1

            # Assign new positions
            newPos[0][i] = pos[0][i] - xcom
            newPos[1][i] = pos[1][i] - ycom
            newPos[2][i] = pos[2][i] - zcom

    # Scaling the cloud to the desired size
    newPos = newPos * (cloudSize/radius)

    # Calculating the volume
    volume = (4. * np.pi /3.) * radius**3

    return ngasNew, newPos, volume

