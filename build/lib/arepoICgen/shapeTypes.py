# Needed libraries
import numpy as np

# Function for a spherical cloud
def sphericalCloud(pos, radius, ngas, verbose=False):
    # Only keeping the particles within the spherical region
    newPos = np.zeros((3, ngas), dtype=np.float64)
    ngasNew = 0

    for i in range(ngas):
        # Distance of each particle from centre
        r = np.sqrt((pos[0,i] - radius)**2 + (pos[1,i] - radius)**2 + (pos[2,i] - radius)**2)

        # Check if within sphere
        if r <= radius:
            # Assign new positions
            newPos[0,ngasNew] = pos[0,i] - radius
            newPos[1,ngasNew] = pos[1,i] - radius
            newPos[2,ngasNew] = pos[2,i] - radius

            # Update the number of particles inside the sphere
            ngasNew += 1

    # Reallocating the position array
    pos = np.zeros((3, ngasNew), dtype=np.float64)
    pos[0] = newPos[0,0:ngasNew]
    pos[1] = newPos[1,0:ngasNew]
    pos[2] = newPos[2,0:ngasNew]

    # Checking the size etc
    xmin = np.min(pos[0])
    xmax = np.max(pos[0])
    ymin = np.min(pos[1])
    ymax = np.max(pos[1])
    zmin = np.min(pos[2])
    zmax = np.max(pos[2])

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

    if verbose:
        # Printing the change in particles
        print("Number originally in box: %s" % ngas)
        print("Number left in sphere: %s" % ngasNew)

        # Printing the new limits
        print("Sphere X Limits: {:.2f} - {:.2f}".format(xmin, xmax))
        print("Sphere Y Limits: {:.2f} - {:.2f}".format(ymin, ymax))
        print("Sphere Z Limits: {:.2f} - {:.2f}".format(zmin, zmax))

        # Printing sphere information
        print("Radius of the sphere: {:.2f}".format(radius))
        print("Sphere COM: {:.2f},{:.2f},{:.2f}".format(xcom, ycom, zcom))
        print("Sphere volume: {:.2e}".format(volume))

    return ngasNew, pos, volume

# Function for an ellipsodial cloud
def ellipsoidalCloud(lengths, ngas, verbose="False"):
    # Generate array for positions
    pos = np.zeros((3, ngas), dtype=np.float64)

    # Calculate radii bounds
    xx = lengths[0]**2
    yy = lengths[1]**2
    zz = lengths[2]**2

    # Generate points on an ellipse
    i = 0
    while i < ngas:
        # Generate the x, y and z coordinates of the points
        x = -lengths[0] + 2 * lengths[0] * np.random.random()
        y = -lengths[1] + 2 * lengths[1] * np.random.random()
        z = -lengths[2] + 2 * lengths[2] * np.random.random()

        # Check we've the right length
        if (x*x/xx + y*y/yy + z*z/zz) <= 1: 
            pos[0,i] = x
            pos[1,i] = y
            pos[2,i] = z
            i += 1

    # Work out volume
    volume = (np.pi * 4/3) * (lengths[0] * lengths[1] * lengths[2])

    return pos, volume

def cylindricalCloud(ngas, radius, lengths, verbose="False"):
    # Setup the positions arrays
    pos = np.zeros((3, ngas), dtype=np.float64)

    # Generate points of the cylnder
    i = 0
    while i < ngas:
        # Generate x (cylinder axis), y and z points of the cylinder
        x = -0.5 * lengths[0] + lengths[0] * np.random.random()
        y = - radius + 2 * radius * np.random.random()
        z = - radius + 2 * radius * np.random.random()

        # Check inside the cylinder
        if np.sqrt(y*y + z*z) <= radius:
            pos[0,i] = x
            pos[1,i] = y
            pos[2,i] = z
            i += 1

    # Work out volume
    volume = np.pi * radius**2 * lengths[0]
    
    return pos, volume