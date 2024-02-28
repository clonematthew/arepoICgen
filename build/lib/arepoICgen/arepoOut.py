# Importing libraries
import numpy as np
from scipy.io import FortranFile
import h5py
import os

# Function to export the created particle data to a usable arepo file
def arepoOut(ngas, pos, vels, pIDs, pMass, pEnergy):
    # Initialising the buffer 
    f = FortranFile("uniformSphere_000", "w")

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
    f.write_record(np.float64(pos.T))

    # Writing out the velocities
    f.write_record(np.float64(vels.T))

    # Writing the particle ids 
    f.write_record(np.int32(pIDs))

    # Writing the particle masses
    f.write_record(np.float64(pMass))

    # Writing the particle internal energies
    f.write_record(np.float64(pEnergy))

    # Writing some random density information
    f.write_record(np.float64(pMass))

# Function to output hdf5 files
def hdf5out(filename, ngas, pos, vels, pIDs, pMass, pEnergy, bField, density=False, pDensity=0):
    # Get path to directory
    dir_path = os.path.dirname(os.path.realpath(__name__))
    
    # Setup file name
    name = dir_path + "/"+ str(filename) + ".hdf5"

    # Opening the ic file
    with h5py.File(name, "w") as icFile:
        # Creating the hdf5 groups
        header = icFile.create_group("Header")
        part0 = icFile.create_group("PartType0")

        # Writing the entries in the header
        numPart = np.array([ngas, 0, 0, 0, 0, 0], dtype=np.int32)
        header.attrs.create("NumPart_ThisFile", numPart)
        header.attrs.create("NumPart_Total", numPart)
        header.attrs.create("NumPart_Total_HighWord", np.zeros(6, dtype=np.int32))
        
        header.attrs.create("MassTable", np.zeros(6, dtype=np.int32))
        header.attrs.create('Time', 0.0)
        header.attrs.create('Redshift', 0.0)
        header.attrs.create('BoxSize', 1.01*np.max(pos))
        header.attrs.create('NumFilesPerSnapshot', 1)
        header.attrs.create('Omega0', 0.0)
        header.attrs.create('OmegaB', 0.0)
        header.attrs.create('OmegaLambda', 0.0)
        header.attrs.create('HubbleParam', 1.0)
        header.attrs.create('Flag_Sfr', 0)
        header.attrs.create('Flag_Cooling', 0)
        header.attrs.create('Flag_StellarAge', 0)
        header.attrs.create('Flag_Metals', 0)
        header.attrs.create('Flag_Feedback', 0)
        header.attrs.create('Flag_DoublePrecision', 1)

        # Extracting the position and velocity components
        x = pos[0]
        y = pos[1]
        z = pos[2]
        vx = vels[0]
        vy = vels[1]
        vz = vels[2]

        # Creating arrays of the right shape to write
        writePos = np.zeros((ngas, 3))
        writeVels = np.zeros((ngas, 3))

        # Assigning the values to the new array in the correct orientation.
        for i in range(ngas):
            writePos[i, 0] = x[i]
            writePos[i, 1] = y[i]
            writePos[i, 2] = z[i]

            writeVels[i, 0] = vx[i]
            writeVels[i, 1] = vy[i]
            writeVels[i, 2] = vz[i]

        # Writing the data of the particles
        part0.create_dataset("ParticleIDs", data=pIDs)
        part0.create_dataset("Coordinates", data=writePos)
        part0.create_dataset("Velocities", data=writeVels)
        part0.create_dataset("InternalEnergy", data=pEnergy)

        # Writing out masses or density based on config
        if density == False:
            # Writing out masses
            part0.create_dataset("Masses", data=pMass)
        elif density == True:
            # Writing out densities
            part0.create_dataset("Masses", data=pDensity)

        # Writing out magnetic field info based on config
        if bField == True:
            # Writing out magnetic field info
            part0.create_dataset("MagneticField", data=np.zeros_like(pMass))
            part0.create_dataset("MagneticFieldDivergence", data=np.zeros_like(pMass))
