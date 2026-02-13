# Importing libraries
import numpy as np
from scipy.io import FortranFile
import h5py
import os

# Function to output hdf5 files
def hdf5out(ngas, pos, vels, pMass, pEnergy, pIDs, config):
    # Get path to directory
    dir_path = os.path.dirname(os.path.realpath(__name__))
    
    # Setup file name
    name = dir_path + "/"+ str(config["filename"]) + ".hdf5"

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
        part0.create_dataset("Masses", data=pMass) # pMass will be density if config["outValue"] = "density"

        # Writing out magnetic field info based on config
        if config["bField"]:
            # Writing out magnetic field info
            part0.create_dataset("MagneticField", data=np.zeros_like(pMass))
            part0.create_dataset("MagneticFieldDivergence", data=np.zeros_like(pMass))
