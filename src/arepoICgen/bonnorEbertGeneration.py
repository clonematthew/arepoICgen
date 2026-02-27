# Library imports
import numpy as np
from random import random
import matplotlib.pyplot as plt

# Function to solve the Lame-Emden equation and work out the density for each cell
def createBonnorEbertSphere(config, params):
    # Calculate the critical radius for the given mass and density
    criticalRadius = (params["mass"]*1.991e33) * 1.66e-24 * 6.67e-8 / (2.4 * 1.38e-16 * params["temp"])
    centralDensity = 18.7 * (1.38e-16**3 * params["temp"]**3) / ((params["mass"]*1.991e33)**2 * 1.66e-24**3 * 6.67e-8**3)
    print("Creating BE sphere with central density {:.2e} gcm/3 and radius {:.0f} AU.".format(centralDensity, criticalRadius/1.5e13))

    # Calculate the radii of each shell and calculate xi
    shellRadii = np.array([0])
    shellRadii = np.append(shellRadii, np.linspace(5e15, criticalRadius, 500))
    beta = 1.38e-16 * params["temp"] / (4 * np.pi * 6.67e-8 * 1.66e-24)
    xi = shellRadii / (beta**0.5 * centralDensity**(-0.5))
    
    # Integrate the Lame-Emden equation to calculate density at every cell point
    psi1 = np.zeros_like(xi)
    psi2 = np.zeros_like(xi)
    psi3 = np.zeros_like(xi)
    dx = np.max(xi) / len(xi)

    for idx, x in enumerate(xi[:-1]):
        psi1[idx+1] = psi1[idx] + dx*psi2[idx]
        psi2[idx+1] = psi2[idx] + dx*ode2(x, psi1[idx], psi2[idx])
        psi3[idx+1] = psi1[idx+1]

    # Calculate density of all cells from the result and scale positions
    density = centralDensity * np.exp(-psi3)

    # Create array for the positions
    positions = np.zeros((3, int(1.3 * params["ngas"])), dtype=np.float64)
    
    # Calculate parameters for populating 
    cellMass = params["mass"] / params["ngas"]
    volume = 4 * np.pi * criticalRadius**3 / 3
    
    # Loop through each of the shells 
    nPlacedCells = 0
    for i in range(len(shellRadii)-1):
        if i == 0:
            # Find the mass in this shell
            shellVolume = 4 * np.pi * shellRadii[i+1]**3 / 3
            shellMass = density[i] * shellVolume
            nCellsInShell = np.round(shellMass / (cellMass * 1.991e33), 0)
            cellSpacing = (shellVolume / nCellsInShell)**(1/3)
            n = int(2 * shellRadii[i+1] / cellSpacing)
            
            # Place the cells in a grid pattern
            from .boxCreation import tripleLoop
            positionsTest = np.zeros((3, n**3))
            positionsTest = tripleLoop(n,n,n, positionsTest, cellSpacing)
            
            for j in range(n**3):
                r = np.sqrt((positionsTest[0,j] - shellRadii[i+1])**2 + (positionsTest[1,j] - shellRadii[i+1])**2 + (positionsTest[2,j] - shellRadii[i+1])**2)
                if r < shellRadii[i+1]:
                    positions[0,nPlacedCells] = positionsTest[0,j] - shellRadii[i+1]/2
                    positions[1,nPlacedCells] = positionsTest[1,j] - shellRadii[i+1]/2
                    positions[2,nPlacedCells] = positionsTest[2,j] - shellRadii[i+1]/2
                    nPlacedCells += 1
        else:    
            # Calculate number of cells we need to populate box with at this density
            nCells = int(np.round(volume * density[i] / (1.991e33*cellMass), 0))

            xs = -criticalRadius + 2. * criticalRadius * np.random.random(nCells)
            ys = -criticalRadius + 2. * criticalRadius * np.random.random(nCells)
            zs = -criticalRadius + 2. * criticalRadius * np.random.random(nCells)
            rs = np.sqrt(xs**2 + ys**2 + zs**2)
                
            # Find cells that are in our shell
            for j in range(int(nCells)):
                if rs[j] > shellRadii[i] and rs[j] < shellRadii[i+1]:
                    positions[0][nPlacedCells] = xs[j]
                    positions[1][nPlacedCells] = ys[j]
                    positions[2][nPlacedCells] = zs[j]
                    nPlacedCells += 1
        
    # Resize position array with amount of cells we placed  
    newPositions = np.zeros((3, nPlacedCells))
    newPositions[0] = positions[0, 0:nPlacedCells]
    newPositions[1] = positions[1, 0:nPlacedCells]
    newPositions[2] = positions[2, 0:nPlacedCells]
                
    return newPositions, nPlacedCells, volume
    
# Lame-Emden function
def ode2(xi, z1, z2):
    if xi == 0:
        return np.exp(-z1)
    else:
        return -((2* z2) / xi) + np.exp(-z1)