######################################
# AREPO Initial Conditions Generator #
# mc 18/10/2023 based on work by pcc #
######################################

# Imports
import numpy as np

# Generate the initial conditions 
def generateICs(config, params):
    # Defining the code units
    uMass = 1.991e33    # grams
    uDist = 1e17        # cm
    uVelo = 36447.2682  # cm/s
    uEner = 1.328e9     # ergs

    # Setting ngas
    ngas = int(params["ngas"])

    #######################
    # Grid type selection #
    #######################

    # Uniform particle grid setups
    if config["grid"] == "boxGrid":
        from .boxCreation import boxGrid

        # Creating a box grid
        pos, ngas, bounds, volume = boxGrid(ngas, params["bounds"])

        # Running the spherical cut module if sphere selected
    elif config["grid"] == "sphereGrid":
        # Import modules for the box and then spherical grid
        from .boxCreation import boxGrid
        from .shapeTypes import sphericalCloud

        # Increase ngas as we will lose some particles when cutting out the sphere
        ngas = int(ngas * 6 / np.pi)

        # Creating a box grid
        pos, ngas, bounds, volume = boxGrid(ngas, params["bounds"])

        # Creating spherical grid
        ngas, pos, volume = sphericalCloud(pos, params["radii"], ngas, bounds)

    # Randomly placed particle setups
    elif config["grid"] == "boxRan":
        from .boxCreation import boxRandom

        # Creating random box grid
        pos, volume = boxRandom(ngas, params["bounds"])

    elif config["grid"] == "sphereRan":
        from .boxCreation import sphereRandom

        # Creating a random spherical grid
        pos, volume = sphereRandom(ngas, params["radii"])

    elif config["grid"] == "ellipseRan":
        from .shapeTypes import ellipsoidalCloud

        # Creating ellipsoid cloud
        pos = ellipsoidalCloud(params["ellipseX"], params["ellipseY"], params["ellipseZ"], ngas)

    # Adjusting positions to be in cm
    pos = pos * 3.09e18

    ###########################
    # Mass and energy defines #
    ###########################

    from .massAndEnergy import masses
    from .massAndEnergy import thermalEnergy

    # Setting equal particle masses
    pMass = masses(ngas, params["mass"])

    # Converting mass into grams
    pMass = pMass * 1.991e33 

    # Working out internal energy of each particle along with the sound speed
    pEnergy, cs = thermalEnergy(ngas, params["temp"], params["mu"])

    # Converting energy into ergs 
    pEnergy = pEnergy * 1e7

    ################################
    # Velocities: Turbulence setup #
    ################################

    # Setup for turbulence from a velocity cube file
    if config["turbulence"] == "turbFile":
        print("Assigning turbulent velocities")
        from .turbulence import turbulenceFromFile

        # Loading in the turbulent velocities from the velocity cube
        velx, vely, velz = turbulenceFromFile(int(config["turbSize"]), config["turbFile"])

        # Branch for the box scenarios
        if config["grid"] == "boxGrid" or config["grid"] == "boxRan":
            from .turbulence import boxGridTurbulence

            # Interpolating and assignning velocities
            vels = boxGridTurbulence(velx, vely, velz, pos, pMass, int(config["turbSize"]), params["virialParam"])

        # Branch for the spherical scenarios
        elif config["grid"] == "sphereGrid" or config["grid"] == "sphereRan":
            from .turbulence import sphericalGridTurbulence

            # Interpolating and assigning velocities
            vels = sphericalGridTurbulence(velx, vely, velz, pos, pMass, int(config["turbSize"]), params["virialParam"])
        else:
            pass
    else:
        # Assgining an empty velocity array if no tubulence setup
        vels = np.zeros((3, ngas), dtype=np.float64)

    ########################
    # Velocities: Rotation #
    ########################

    # Add rotation to the body
    if config["rotation"] == "rotation":
        print("Adding solid body rotation")
        from .rotation import addRotation

        # Add rotation around z axis of given beta energy ratio
        vels = addRotation(pos, pMass, vels, params["beta"])
    else:
        pass

    #####################
    # Special Functions #
    #####################

    # Add a Boss-Bodenheimer density perturbation (Boss & Bodenheimer 1979)
    if config["extras"] == "bossBodenheimer":
        print("Adding Boss-Bodenheimer perturbation")
        from .densityPerturbations import bossBodenheimer
        pos, pMass = bossBodenheimer(ngas, pos, pMass)

    ###################################
    # Setting particle identification #
    ###################################

    # Assigning each particle an ID from 1 to the max number of particles
    pIDs = np.linspace(1, ngas, ngas, dtype=np.int32)

    ################################
    # Low density particle padding #
    ################################

    # Setup for padding the box with low density particles
    if config["padding"] == True:
        # Branch for the box setups
        print("Padding box with low density particles")
        if config["grid"] == "boxGrid" or config["grid"] == "boxRan":
            from .lowDensityPadding import padBox

            # Pad the box with low density particles outside the box grid
            pos, vels, pMass, pIDs, pEnergy, pRho, ngasAll = padBox(ngas, pos, vels, pMass, pIDs, pEnergy, params["boxDims"], params["tempFactor"])
        
        # Branch for the spherical setups
        elif config["grid"] == "sphereGrid" or config["grid"] == "sphereRan":
            from .lowDensityPadding import padSphere

            # Pad the box with low density particles outside the spherical cloud
            pos, vels, pMass, pIDs, pEnergy, pRho, ngasAll = padSphere(ngas, pos, vels, pMass, pIDs, pEnergy, params["boxDims"], params["tempFactor"])
        else:
            ngasAll = ngas
    else:
        ngasAll = ngas
        pass

    #############################
    # Making positions positive #
    #############################

    # Getting the minimum value of every coordinate
    minx = np.min(pos[0])
    miny = np.min(pos[1])
    minz = np.min(pos[2])

    # Shifting everything if its less than zero
    if minx < 0:
        pos[0] -= minx
    if miny < 0:
        pos[1] -= miny
    if minz < 0:
        pos[2] -= minz

    ############################################
    # Conversion of quantities into code units #
    ############################################

    # All variables should be in c.g.s units for conversion
    pos = pos / uDist
    vels = vels / uVelo
    pMass = pMass / uMass
    pEnergy = pEnergy / uEner

    ##############################
    # Desired Density Conversion #
    ##############################

    if config["outValue"] == "density":
        # Converting the number density to code units
        densityTarget = params["density"] * params["mu"] * 1.66e-24 
        densityTarget = densityTarget / (uMass / (uDist**3))
        densityTargetPadding = densityTarget * 0.01

        # Creating density array
        pDensity = np.ones_like(pMass)
        pDensity[0:ngas] = densityTarget
        pDensity[ngas:-1] = densityTargetPadding

    ########################
    # File output to AREPO #
    ########################

    print("Writing output file")

    if config["output"] == "hdf5":
        from .arepoOut import hdf5out

        # Writing masses to mass 
        if config["outValue"] == "masses":
            # Write the particle data as a hdf5 file
            hdf5out(config["filename"], ngasAll, pos, vels, pIDs, pMass, pEnergy, config["bField"])

        # Writing density to mass
        elif config["outValue"] == "density":
            # Write the particle data as a hdf5 file
            hdf5out(config["filename"], ngasAll, pos, vels, pIDs, pMass, pEnergy, config["bField"], True, pDensity)
    else:
        print("Fortran binary version is broken, sorry </3")

# Function to just easily create a uniform sphere       
def easySphere(mass, numDense, ngas, filename, mu=1.4, beta=0.05):
    # Convert number density to density
    density = numDense * 1.66e-24 * mu

    # Calculate the radius of the sphere
    volume = (mass*1.991e33) / density
    radius = (volume * 3/(4*np.pi))**(1/3)
    radius = radius / 3.09e18

    # Generate the initial conditions
    config = {
        "grid": "sphereGrid",
        "turbulence": "static",
        "rotation": "rotation",
        "padding": True,
        "output": "hdf5",
        "outValue": "masses",
        "extras": "none",
        "bField": False,
        "filename": filename
    }

    params = {
        "ngas": ngas,
        "bounds": [0, radius*5, 0, radius*5, 0, radius*5],
        "radii": radius,
        "mass": mass,
        "temp": 15,
        "mu": 1.4,
        "beta": beta,
        "boxDims": [5, 5, 5],
        "tempFactor": 2,
    }

    generateICs(config, params)