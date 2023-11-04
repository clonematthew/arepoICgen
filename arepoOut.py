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
    nPartHW = np.zeros(6, dtype=np.int32)

    # Setting the mass array
    massArray = np.zeros(6, dtype=np.float64)

    # Setting other variables
    time = np.float64(0.)
    redshift = np.float64(0.)
    sfrFlag = np.int32(0)
    feedbackFlag = np.int32(0)

    nAll = np.zeros(6, dtype=np.int32)
    nAll[0] = np.int32(ngas)

    coolingFlag = np.int32(0)

    numFiles = np.int32(1)
    boxsize = np.float64(0.)
    omega0 = np.float64(0.)
    omegaLambda = np.float64(0.)
    hubbleParam = np.float64(0.)
    stellarangeFlag = np.int32(0)
    metalsFlag = np.int32(0)
    entropyicFlag = np.int32(0)

    doublePrescisionFlag = np.int32(1)
    lptICsFlag = np.int32(0)
    lptScalingFactor = np.int32(0)
    tracerFieldFlag = np.int32(0)
    compositionVectorLength = np.int32(0)

    # Defining the variable to fill up the rest of the buffer
    unused = np.zeros(10, dtype=np.int32)

    # Writing all of the header information
    f.write_record(nPartArray, massArray, time, redshift, sfrFlag, feedbackFlag, nAll, coolingFlag, 
                   numFiles, boxsize, omega0, omegaLambda, hubbleParam, stellarangeFlag, metalsFlag, 
                   nPartHW, entropyicFlag, doublePrescisionFlag, lptICsFlag, lptScalingFactor, 
                   tracerFieldFlag, compositionVectorLength, unused)
    
    # Writing out the positions
    f.write_record(np.float64(pos))

    # Writing out the velocities
    f.write_record(np.float64(vels))

    # Writing the particle ids 
    f.write_record(np.int32(pIDs))

    # Writing the particle masses
    f.write_record(np.float64(pMass))

    # Writing the particle internal energies
    f.write_record(np.float64(pEnergy))


