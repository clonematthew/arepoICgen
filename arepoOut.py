# Importing libraries
import numpy as np
from scipy.io import FortranFile

# Function to export the created particle data to a usable arepo file
def arepoOut(ngas, pos, vels, pIDs, pMass, pEnergy):
    # Initialising the buffer 
    f = FortranFile("arepo_out", "w")

    # Setting the number of particles
    nPartArray = np.zeros(6, dtype=np.int32)
    nPartArray[0] = np.int32(ngas)

    # Setting the mass array
    massArray = np.zeros(6, dtype=np.double)

    # Setting other variables
    time = np.double(0.)
    redshift = np.double(0.)
    sfrFlag = np.int32(0)
    feedbackFlag = np.int32(0)

    nAll = np.zeros(6, dtype=np.int32)
    nAll[0] = np.int32(ngas)

    coolingFlag = np.int32(0)

    numFiles = np.int32(1)
    boxsize = np.double(0.)
    omega0 = np.double(0.)
    omegaLambda = np.double(0.)
    hubbleParam = np.double(0.)
    stellarangeFlag = np.int32(0)
    metalsFlag = np.int32(0)
    entropyicFlag = np.int32(0)

    # Defining the variable to fill up the rest of the buffer
    unused = np.zeros(12, dtype=np.int32)

    # Writing all of the header information
    f.write_record(nPartArray, massArray, time, redshift, sfrFlag, feedbackFlag, nAll, coolingFlag, numFiles, 
                  boxsize, omega0, omegaLambda, hubbleParam, stellarangeFlag, metalsFlag, nAll,
                  np.int32(0), np.int32(0), np.float32(0.), entropyicFlag, unused)
    
    # Writing out the positions
    f.write_record(np.double(pos))

    # Writing out the velocities
    f.write_record(np.double(vels))

    # Writing the particle ids 
    f.write_record(np.int32(pIDs))

    # Writing the particle masses
    f.write_record(np.double(pMass))

    # Writing the particle internal energies
    f.write_record(np.double(pEnergy))


