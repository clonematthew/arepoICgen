# Importing libraries
import numpy as np
from scipy.io import FortranFile

# Function to export the created particle data to a usable arepo file
def arepoOut(ngas, pos, vels, pIDs, pMass, pEnergy):
    # Starting with the header file, defining all the characteristics 
    nPartArray = np.zeros(6, dtype=np.int32)
    nAll = np.zeros(6, dtype=np.int32)

    # Setting the number of particles
    nPartArray[0] = ngas
    nAll[0] = ngas
    nTotal = np.sum(nAll, dtype=np.int32)

    # Setting the mass array
    massArray = np.double(0.)

    # Setting other variables
    time = np.double(0.)
    redshift = np.double(0.)

    sfrFlag = np.int32(0)
    feedbackFlag = np.int32(0)
    coolingFlag = np.int32(0)
    stellarangeFlag = np.int32(0)
    metalsFlag = np.int32(0)
    entropyicFlag = np.int32(0)

    numFiles = np.int32(1)
    boxsize = np.double(0.)
    omega0 = np.double(0.)
    omegaLambda = np.double(0.)
    hubbleParam = np.double(0.)

    # Unit variables
    uMass = np.float64(1.991e33)
    uDist = np.float64(1.000e17)
    uTime = np.float64(0.)

    # Defining the variable to fill up the rest of the buffer
    unused = np.zeros(24, dtype=np.int32)

    # Initialising the buffer 
    f = FortranFile("arepo_out", "w")

    # Writing all of the header information
    f.write_record(nPartArray, massArray, time, redshift, sfrFlag, feedbackFlag, nAll, coolingFlag, numFiles, 
                   boxsize, omega0, omegaLambda, hubbleParam, stellarangeFlag, metalsFlag, nTotal, uTime, uMass,
                    uDist, entropyicFlag, unused)
    
    # Writing out the positions
    f.write_record(pos)

    # Writing out the velocities
    f.write_record(vels)

    # Writing the particle ids 
    f.write_record(pIDs)

    # Writing the particle masses
    f.write_record(pMass)

    # Writing the particle internal energies
    f.write_record(pEnergy)


