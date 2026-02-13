# Needed libraries
import numpy as np

def generateCellGrid(config, params):
    if config["grid"] == "boxGrid":
        from .boxCreation import boxGrid

        # Creating a box grid
        pos, ngas, volume = boxGrid(int(params["ngas"]), params["lengths"], config["verbose"])

        # Running the spherical cut module if sphere selected
    elif config["grid"] == "sphereGrid":
        # Import modules for the box and then spherical grid
        from .boxCreation import boxGrid
        from .shapeTypes import sphericalCloud

        # Increase ngas as we will lose some particles when cutting out the sphere
        ngas = int(int(params["ngas"]) * 6 / np.pi)

        # Our box will always be 2x the radius in each dimension
        params["lengths"] = [2*params["radii"], 2*params["radii"], 2*params["radii"]]

        # Creating a box grid
        pos, ngas, volume = boxGrid(ngas, params["lengths"], config["verbose"])

        # Cutting a sphere out of it 
        ngas, pos, volume = sphericalCloud(pos, params["radii"], ngas, config["verbose"])

    # Randomly placed particle setups
    elif config["grid"] == "boxRan":
        from .boxCreation import boxRandom

        # Creating random box grid
        ngas, pos, volume = boxRandom(int(params["ngas"]), params["lengths"], config["verbose"])

    elif config["grid"] == "sphereRan":
        from .boxCreation import sphereRandom

        # Creating a random spherical grid
        ngas, pos, volume = sphereRandom(int(params["ngas"]), params["radii"], config["verbose"])

    elif config["grid"] == "ellipseRan":
        from .shapeTypes import ellipsoidalCloud

        # Creating ellipsoid cloud
        ngas, pos, volume = ellipsoidalCloud(params["lengths"], int(params["ngas"]), config["verbose"])

    elif config["grid"] == "cylinderRan":
        from .shapeTypes import cylindricalCloud

        # Creating cylinderical cloud
        ngas, pos, volume = cylindricalCloud(int(params["ngas"]), params["radii"], params["lengths"], config["verbose"])
    
    # Ensure that all the setups are centered on the 0,0 axis
    if np.min(pos[0]) >= 0:
        pos[0] -= (np.max(pos[0]) - np.min(pos[0]))
    if np.min(pos[1]) >= 0:
        pos[1] -= (np.max(pos[1]) - np.min(pos[1]))
    if np.min(pos[2]) >= 0:
        pos[2] -= (np.max(pos[2]) - np.min(pos[2]))
    
    return ngas, pos, volume