'''
                                                  /$$$$$$  /$$$$$$                               
                                                 |_  $$_/ /$$__  $$                              
  /$$$$$$   /$$$$$$   /$$$$$$   /$$$$$$   /$$$$$$  | $$  | $$  \__/  /$$$$$$   /$$$$$$  /$$$$$$$ 
 |____  $$ /$$__  $$ /$$__  $$ /$$__  $$ /$$__  $$ | $$  | $$       /$$__  $$ /$$__  $$| $$__  $$
  /$$$$$$$| $$  \__/| $$$$$$$$| $$  \ $$| $$  \ $$ | $$  | $$      | $$  \ $$| $$$$$$$$| $$  \ $$
 /$$__  $$| $$      | $$_____/| $$  | $$| $$  | $$ | $$  | $$    $$| $$  | $$| $$_____/| $$  | $$
|  $$$$$$$| $$      |  $$$$$$$| $$$$$$$/|  $$$$$$//$$$$$$|  $$$$$$/|  $$$$$$$|  $$$$$$$| $$  | $$
 \_______/|__/       \_______/| $$____/  \______/|______/ \______/  \____  $$ \_______/|__/  |__/
                              | $$                                  /$$  \ $$                    
                              | $$                                 |  $$$$$$/                    
                              |__/                                  \______/                     
               
    arepoICgen: Function to create initial conditions for the moving mesh code AREPO.  
                Written by Matt Cusack, 2025. Based on setpartarepo by Paul Clark.             
'''

# Library imports
import numpy as np

# Main class to generate ICs
class arepoICgen():
    # Initialise the class and check the input settings
    def __init__(self, config, params):
        # Validate the input settings
        self.validateInput(config, params)
        
        # Generate the initial conditions using given settings
        self.generateICs()
        
        # Output the initial conditions to a hdf5 file
        self.outputICs()
        
    # Main function to generate the initial conditions
    def generateICs(self):
        # Generate the cell grid
        from .gridGenerator import generateCellGrid
        self.ngas, self.positions, self.volume = generateCellGrid(self.config, self.params)
        
        # Set mass and temperature of the cells
        from .massAndEnergy import masses, thermalEnergy
        self.masses = masses(self.ngas, self.config, self.params) * 1.991e33
        self.energies = thermalEnergy(self.ngas, self.config, self.params) * 1e7
        
        # Add velocity to the cells
        from .velocityGenerator import addVelocities
        self.velocities = addVelocities(self.positions, self.masses, self.energies, self.config, self.params)
        
        # Make modifications to the ICs
        if self.config["extras"] != "none":
            self.modifyICs()
        
        # Pad the box around the cloud
        from .boxPadding import padBox
        self.ngas, self.positions, self.velocities, self.masses, self.energies = padBox(self.ngas, self.positions, self.velocities,
                                                                                        self.masses, self.energies, self.volume,
                                                                                        self.config, self.params)
            
        # Provide cells with unique IDs
        self.cellIDs = np.linspace(1, self.ngas, self.ngas, dtype=np.int32)
    
    # Function to modify the ICs using different density perturbations etc
    def modifyICs(self):
        # Add a Boss-Bodenheimer density perturbation (Boss & Bodenheimer 1979)
        if self.config["extras"] == "bossBodenheimer":
            from .densityPerturbations import bossBodenheimer
            self.positions, self.masses = bossBodenheimer(self.ngas, self.positions, self.masses)
            
        # Add a density gradient along one axis
        elif self.config["extras"] == "densityGradient":
            from .densityPerturbations import densityGradient
            self.masses = densityGradient(self.positions, self.masses, self.params)
            
        # Add a Bonnor-Ebert density profile (Bonnor 1956, Ebert 1955)
        elif self.config["extras"] == "bonnorEbert":
            from .bonnorEbertGeneration import createBonnorEbertSphere
            self.positions, self.masses = createBonnorEbertSphere(self.ngas, self.positions, self.params)
        
    # Function to output the initial conditions to a file
    def outputICs(self):
        # Shift all cells such that they are all > 0
        if np.min(self.positions[0]) < 0:
            self.positions[0] -= np.min(self.positions[0])
        if np.min(self.positions[1]) < 0:
            self.positions[1] -= np.min(self.positions[1])
        if np.min(self.positions[2]) < 0:
            self.positions[2] -= np.min(self.positions[2])
            
        # Convert all variables into code units
        self.positions = self.positions / self.uDist
        self.velocities = self.velocities / self.uVelo
        self.energies = self.energies / self.uEner
        
        # If we're outputting density as mass, convert accordingly
        if self.config["outValue"] == "density":
            self.masses = self.masses / (self.uMass / self.uDist**3)
        else:
            self.masses = self.masses / self.uMass       
            
        # Output the data to a hdf5 file
        from .arepoOut import hdf5out
        hdf5out(self.ngas, self.positions, self.velocities, self.masses, self.energies, self.cellIDs, self.config)   
        
    # Function to pass through configs and params and check if valid/present
    def validateInput(self, config, params):
        # Assign input settings to class
        self.config = config
        self.params = params
        
        # Get the keys of each of the settings dictionaries
        allConfigs = config.keys()
        allParams = params.keys()
        
        # Check grid settings
        if "ngas" not in allParams:
            raise Exception("No ngas! How many cells do you want?")
            
        if "grid" not in allConfigs:
            raise Exception("No grid! How do you want the cells laid out?")
        else:
            if config["grid"] == "sphereRan" or config["grid"] == "sphereGrid":
                if "radii" not in allParams:
                    raise Exception("No radius! How big is the sphere supposed to be?")
            else:
                if "lengths" not in allParams:
                    raise Exception("No box side lengths! How big is the box?")
                    
        # Check mass and energy settings
        if "mass" not in allParams:
            raise Exception("No mass! How massive is the cloud?")
        if "temp" not in allParams:
            print("No temperature! Assuming T = 20K.")
            params["temp"] = 20
        if "mu" not in allParams:
            print("No mu. Assuming mu = 1.4.")
            params["mu"] = 1.4
            
        # Check velocity settings
        if "turbulence" not in allConfigs:
            print("No turbulence specified, using static cloud.")
            config["turbulence"] = "static"
        elif config["turbulence"] == "turbFile":
            if "turbFile" not in allConfigs:
                raise Exception("No turbulence file! No grid to interpolate from.")
            if "turbSize" not in allConfigs:
                print("No turb file size. Assuming 128x128x128 grid.")
                config["turbSize"] = 128
            if "epsilon" not in allParams:
                print("No virial parameter specified. Assuming virialised cloud.")
                params["virialParam"] = 2
        if "rotation" not in allConfigs:
            config["rotation"] = "none"
        elif config["rotation"] == "rotation":
            if "beta" not in allParams:
                print("No beta specified. Assuming 0.01.")
                params["beta"] = 0.01
                
        # Check padding settings
        if "padBox" not in allConfigs:
            config["padBox"] = False
        else:
            if "paddingPercent" not in allParams:
                print("No padding particles percent specified, using 0.02.")
                params["paddingPercent"] = 0.02
            if "paddingDensity" not in allParams:
                print("No padding density specified, using 0.01.")
                params["paddingDensity"] = 0.01
            if "boxSize" not in allParams:
                print("No padding box dimensions provided, using 2x.")
                params["boxSize"] = [2, 2, 2]
            if "tempFactor" not in allParams:
                print("No temperature factor for padding cells provided, using 2x.")
                params["tempFactor"] = 2      
                
        # Check output settings
        if "outValue" not in allConfigs:
            config["outValue"] = "mass"
        if "filename" not in allConfigs:
            raise Exception("No filename given. Please provide filename to save the intitial conditions to.")
        if "bField" not in allConfigs:
            config["bField"] = False
            
        # Check extras settings                    
        if "extras" not in allConfigs:
            config["extras"] = "none"
        else:
            if config["extras"] == "bonnorEbert" and config["outValue"] != "density":
                print("Bonnor-Ebert sphere needs to output as density, forcing density output.")
                config["outValue"] = "density"
            
        # Check verbose setting
        if "verbose" not in allConfigs:
            config["verbose"] = False
                                
        # Set default code units
        if "uMass" not in allParams:
            self.uMass = 1.991e33
        else:
            self.uMass = params["uMass"]
        if "uDist" not in allParams:
            self.uDist = 1e17 
        else:
            self.uDist = params["uDist"]
        if "uVelo" not in allParams:
            self.uVelo = 36447.2682
        else:
            self.uVelo = params["uVelo"]
        if "uEner" not in allParams:
            self.uEner = 1.328e9
        else:
            self.uEner = params["uEner"]  