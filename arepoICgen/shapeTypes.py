# Needed libraries
import numpy as np

# Function for a spherical cloud
def sphericalCloud(pos, cloudSize, ngas, bounds):
    # Getting the min and max positions
    xmin = bounds[0]
    xmax = bounds[1]
    ymin = bounds[2]
    ymax = bounds[3]
    zmin = bounds[4]
    zmax = bounds[5]

    # Getting box sizes
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin

    # Finding the radius of the central sphere
    radius = np.abs(np.min([dx/2, dy/2, dz/2]))

    # Fiding the central point
    xcom = xmin + dx/2
    ycom = ymin + dy/2
    zcom = zmin + dz/2

    # Only keeping the particles within the spherical region
    newPos = np.zeros((3, ngas), dtype=np.float64)
    ngasNew = 0

    for i in range(ngas):
        # Distance of each particle from centre
        r = np.sqrt((pos[0,i] - xcom)**2 + (pos[1,i] - ycom)**2 + (pos[2,i] - zcom)**2)

        # Check if within sphere
        if r <= radius:
            # Assign new positions
            newPos[0,ngasNew] = pos[0,i] - xcom
            newPos[1,ngasNew] = pos[1,i] - ycom
            newPos[2,ngasNew] = pos[2,i] - zcom

            # Update the number of particles inside the sphere
            ngasNew += 1

    # Printing the change in particles
    print("Number originally in box: %s" % ngas)
    print("Number left in sphere: %s" % ngasNew)

    # Reallocating the position array
    pos = np.zeros((3, ngasNew), dtype=np.float64)
    pos[0] = newPos[0,0:ngasNew]
    pos[1] = newPos[1,0:ngasNew]
    pos[2] = newPos[2,0:ngasNew]

    # Scaling the cloud to the desired size
    pos = pos * (cloudSize/radius)

    # Checking the size etc
    xmin = np.min(pos[0])
    xmax = np.max(pos[0])
    ymin = np.min(pos[1])
    ymax = np.max(pos[1])
    zmin = np.min(pos[2])
    zmax = np.max(pos[2])

    # Printing the new limits
    print("Sphere X Limits: {:.2f} - {:.2f}".format(xmin, xmax))
    print("Sphere Y Limits: {:.2f} - {:.2f}".format(ymin, ymax))
    print("Sphere Z Limits: {:.2f} - {:.2f}".format(zmin, zmax))

    # Getting box sizes
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin

    # Finding the radius of the central sphere
    radius = np.abs(np.min([dx/2, dy/2, dz/2]))

    # Fiding the central point
    xcom = xmin + dx/2
    ycom = ymin + dy/2
    zcom = zmin + dz/2
    
    # Calculating the volume
    volume = (4. * np.pi /3.) * radius**3

    # Printing sphere information
    print("Radius of the sphere: {:.2f}".format(radius))
    print("Sphere COM: {:.2f},{:.2f},{:.2f}".format(xcom, ycom, zcom))
    print("Sphere volume: {:.2e}".format(volume))

    return ngasNew, pos, volume

