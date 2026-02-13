# Library imports
import numpy as np
import matplotlib.pyplot as plt

# Function to solve the Lame-Emden equation and work out the density for each cell
def createBonnorEbertSphere(ngas, positions, params):
    # Calculate the critical radius for the given mass and density
    criticalRadius = 1.8 * (1.38e-16 * params["temp"] / (1.66e-24 * 6.67e-8 * params["centralDensity"]))**0.5

    # Calculate radius of every cell to the centre
    radius = np.sqrt(positions[0]**2 + positions[1]**2 + positions[2]**2) 
    
    # Scale the radius so that its max is the critical radius
    scalingFactor = criticalRadius / np.max(radius)
    radius *= scalingFactor
    radiusSorted = np.argsort(radius)
    xi = radius[radiusSorted] / (criticalRadius / 1.8)
    
    # Integrate the Lame-Emden equation to calculate density at every cell point
    psi1 = np.zeros_like(radius)
    psi2 = np.zeros_like(radius)
    psi3 = np.zeros_like(radius)
    dx = 1.8 / ngas

    for idx, x in enumerate(xi[:-1]):
        psi1[idx+1] = psi1[idx] + dx*psi2[idx]
        psi2[idx+1] = psi2[idx] + dx*ode2(x, psi1[idx], psi2[idx])
        psi3[radiusSorted[idx+1]] = psi1[idx+1]

    # Calculate density of all cells from the result and scale positions
    density = params["centralDensity"] * np.exp(-psi3)
    positions *= scalingFactor
    
    return positions, density
    
# Lame-Emden function
def ode2(xi, z1, z2):
    if xi == 0:
        return np.exp(-z1)
    else:
        return -((2* z2) / xi) + np.exp(-z1)