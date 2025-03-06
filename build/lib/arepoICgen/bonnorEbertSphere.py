# Import libraries
import numpy as np
from numba import jit

# Some useful constants
G = 6.67e-8
kB = 1.38e-16

# Solve the LE equation for an unstable sphere
def solveLEequation(dimensionlessRadius):
    # Work out number of integration steps to do
    dxIntegration = 0.00001
    numberOfPoints = int((dimensionlessRadius/dxIntegration) + 1.1)
    dxHalfStep = dxIntegration / 2

    # Create arrays 
    massFunc = np.zeros(numberOfPoints)
    pressureFunc = np.zeros(numberOfPoints)
    densityFunc = np.zeros(numberOfPoints)
    boundaryRadiusFunc = np.zeros(numberOfPoints)
    radiusFunc = np.zeros(numberOfPoints)

    # Set input values
    radiusFunc[1] = dxIntegration
    massFunc[1] = (dxIntegration**3/3.) - (dxIntegration**5/30.) + (dxIntegration**7/315.)
    densityFunc[1] = (dxIntegration**2/6.) - (dxIntegration**4/120.) + (dxIntegration**6/1890.)

    # Solve the equation using an RK4 algorithim
    for i in range(2, numberOfPoints):
        xi1=radiusFunc[i-1]                                             
        mu1=massFunc[i-1]                                             
        ps1=densityFunc[i-1]                                             
        dmu1=xi1**2*np.exp(-ps1)                                 
        dps1=mu1/(xi1**2)                                    

        xi2=xi1+dxHalfStep                                       
        mu2=mu1+(dmu1*dxHalfStep)                                  
        ps2=ps1+(dps1*dxHalfStep)                                  
        dmu2=xi2**2*np.exp(-ps2)                                
        dps2=mu2/(xi2**2)                                      

        xi3=xi2                                                
        mu3=mu1+(dmu2*dxHalfStep)                                  
        ps3=ps1+(dps2*dxHalfStep)                                 
        dmu3=xi3**2*np.exp(-ps3)                                   
        dps3=mu3/(xi3**2)                                       

        xi4=xi2+dxHalfStep                                     
        mu4=mu1+(dmu3*dxIntegration)                                     
        ps4=ps1+(dps3*dxIntegration)                                     
        dmu4=xi4**2*np.exp(-ps4)                                   
        dps4=mu4/(xi4**2)                                      

        # Assign computed final values for the step
        radiusFunc[i] = xi4                                              
        massFunc[i] = massFunc[i-1]+(dmu1+(2.*dmu2)+(2.*dmu3)+dmu4)*dxIntegration/6.   
        densityFunc[i] = densityFunc[i-1]+(dps1+(2.*dps2)+(2.*dps3)+dps4)*dxIntegration/6.   
        pressureFunc[i] = massFunc[i]**2*np.exp(-densityFunc[i])                             
        boundaryRadiusFunc[i] = radiusFunc[i]/massFunc[i]  
        
        # The edge of the sphere is where the pressure begins to drop
        if pressureFunc[i] < pressureFunc[i-1]:
            # Store critical value
            iCRIT = i-1
            
            # Density contrast at the boundary, ln(rhoCent / rhoBoundary)
            densityContrast = densityFunc[iCRIT]
            
            # Dimenionless mass of the sphere
            dimensionlessMass = massFunc[iCRIT]
            
            aux0 = radiusFunc[iCRIT]**3 * np.exp(-densityContrast) / (3.*dimensionlessMass)
            
    # Create lookup table arrays
    interpolatedMass = 0.001 * np.arange(-9500, 0, 1)
    interpolatedRadius = np.zeros(9500)
    
    massFunc *= 2
    
    # Initialise x and y old values to very low
    xOLD = -99
    yOLD = -99
    
    # Work out auxilliary variables
    aux0 = 1./massFunc[iCRIT]
    aux1 = 1./radiusFunc[iCRIT]
    j = 0
    
    # Loop over the integrated values
    for i in range(1, numberOfPoints):
        xNEW = np.log10(aux0 * massFunc[i])
        yNEW = np.log10(aux1 * radiusFunc[i])
        
        # Interpolate required radius value
        if xNEW > interpolatedMass[j]:
            interpolatedRadius[j] = yOLD + ((interpolatedMass[j] - xOLD) * (yNEW - yOLD) / (xNEW - xOLD))
            
            j += 1
        
        # Update x and y new
        xOLD = xNEW
        yOLD = yNEW
        
        if j == 9500:
            break  
        
    interpolatedRadius[-1] = 0
        
    return radiusFunc, massFunc, densityFunc, densityContrast, dimensionlessMass, interpolatedRadius

# Loop to place particles with (for speed)
@jit(nopython=True)
def particleLoop(nn0, rBoundary, halfWidthOfDomain, rotMatrix, Cmm, interpolatedRadius, ffB, fffB, tempFactor):
    # Setup positions of the cube
    spacing = (np.sqrt(2)/nn0)**0.3333333

    dxi0 = spacing
    iMID = int(halfWidthOfDomain/dxi0) + 2
    iMAX = 2 * iMID - 1

    dxi = np.zeros(iMAX)

    for i in range(iMAX):
        dxi[i] = (i+1) * dxi0
        
    dxj0 = spacing/2
    dyj0 = np.sqrt(3)*spacing/2
    jMID = int(halfWidthOfDomain/dyj0)+2
    jMAX = 2 * jMID - 1

    dxj = np.zeros(jMAX)
    dyj = np.zeros(jMAX)

    for j in range(jMAX):
        dxj[j] = ((j+1 - jMID) % 2) * dxj0
        dyj[j] = (j+1) * dyj0
        
    dxk0 = spacing/4
    dyk0 = spacing/np.sqrt(12)
    dzk0 = spacing/np.sqrt(1.5)
    kMID = int(halfWidthOfDomain/dzk0) + 2
    kMAX = 2 * kMID - 1

    dxk = np.zeros(kMAX)
    dyk = np.zeros(kMAX)
    dzk = np.zeros(kMAX)

    for k in range(kMAX):
        i = 1 - (-1)**((k+1 - kMID) % 3)
        dxk[k] = i * dxk0
        dyk[k] = ((k+1 - kMID) % 3) * dyk0
        dzk[k] = (k+1) * dzk0

    # Calculate offsets    
    dxS = - dxi[iMID] - dxj[jMID] - dxk[kMID]
    dyS = - dyj[jMID] - dyk[kMID]
    dzS = - dzk[kMID]
    
    # Calculate total number of particles
    pTotal = iMAX * jMAX * kMAX

    # Create array of positions and temperatures
    r = np.zeros((3,pTotal))
    t = np.zeros(pTotal)
    
    p = 0
    pInSphere = 0

    rB3 = rBoundary**3
    X03 = halfWidthOfDomain**3
    
    # Loop through each axis
    for i in range(iMAX):
        for j in range(jMAX):
            for k in range(kMAX):
                # Compute x position
                r[0,p] = dxi[i] + dxj[j] + dxk[k] + dxS
                
                # Check if outside domain
                if (abs(r[0,p]) < halfWidthOfDomain):
                    # Compute y position
                    r[1,p] = dyj[j] + dyk[k] + dyS
                    
                    # Check if outside domain
                    if (abs(r[1,p]) < halfWidthOfDomain):
                        # Compute z position
                        r[2,p] = dzk[k] + dzS
                        
                        # Check if outisde domain
                        if (abs(r[2,p]) < halfWidthOfDomain):
                            # Rotate the particle by the random matrix  
                            r[:,p] = rotMatrix[0]*r[0,p] + rotMatrix[1]*r[1,p] + rotMatrix[2]*r[2,p]
                            
                            # Find the radius of this particle
                            rrMAG = r[0,p]**2 + r[1,p]**2 + r[2,p]**2
                            
                            if rrMAG > 0.1e-10:
                                rrMAG = np.sqrt(rrMAG)
                            
                                # Calculate the enclosed mass
                                mEnc = Cmm * rrMAG**3
                            
                                # Inside the BE sphere
                                if mEnc < 1:
                                    pInSphere += 1
                                    
                                    # Calculate T index for lookup
                                    T = 1000 * np.log10(mEnc)
                                    iLO = int(np.floor(T))
                                    iHI = int(np.ceil(T))
                                    
                                    # Calculate radius post-stretch
                                    rMAG = rBoundary * 10**((T - iLO) * interpolatedRadius[iHI+9499] + (iHI - T) * interpolatedRadius[iLO+9499])
                                    t[p] = 1
                                else:
                                    rMAG = (rB3 + fffB * ((ffB * rrMAG**3) - X03))**0.3333333
                                    t[p] = tempFactor
                                
                                # Scale the radius position
                                r[0,p] = (rMAG / rrMAG) * r[0,p]
                                r[1,p] = (rMAG / rrMAG) * r[1,p]
                                r[2,p] = (rMAG / rrMAG) * r[2,p]
                                
                            if abs(r[0,p]) < halfWidthOfDomain and abs(r[1,p]) < halfWidthOfDomain and abs(r[2,p]) < halfWidthOfDomain:
                                p += 1  
                                
    return r, p, pInSphere, t

# Create a BE sphere
def createBEsphere(BEmass, ngas, temperature, mu, paddingDensityContrast, tempFactor, dimensionlessRadius=6.5):
    # Solve the LE equation
    radiusFunc, massFunc, densityFunc, densityContrast, dimensionlessMass, interpolatedRadius = solveLEequation(dimensionlessRadius)
    
    # Calculate sound speed of the sphere
    cs = np.sqrt(kB * temperature / (mu * 1.66e-24))
    
    # Calculate the boundary radius in pc
    rBoundary = 0.0043016 * (BEmass/2) * dimensionlessRadius / ((cs/1e5)**2 * dimensionlessMass)  
    
    # Calculate the number of particles to pad the box with
    ppLag = dimensionlessRadius**3 * np.exp(-densityContrast) * ngas / (12.566371 * dimensionlessMass * (1/paddingDensityContrast))
    fB = 1 + (2 * 10 * 1.2 / (ppLag**0.33333))
    halfWidthOfDomain = fB * rBoundary
    nPaddingParticles = int(((8 * fB**3) - 4.1887902) * ppLag)

    # Estimate the total number of particles
    nTotal = ngas + nPaddingParticles

    # Work out analytical number of particles
    ffB = 1 + (((5.1961524 * fB**3 - 1.) * dimensionlessRadius**3 * np.exp(-densityContrast)) / (3 * dimensionlessMass * (1/paddingDensityContrast)))
    ppTOT = int(1.9098593 * ngas * ffB)

    # Work out number density of particles
    nn0 = 0.23873241 * ffB * ngas / (halfWidthOfDomain**3)

    # Work out pre-stretch radius of the sphere
    rrB = halfWidthOfDomain / (ffB**0.33333)

    # Work out mass coefficient inside radius
    Cmm = ffB / (halfWidthOfDomain**3)
    fffB = (3 * (1/paddingDensityContrast) * dimensionlessMass) / ((fB * dimensionlessRadius)**3 * np.exp(-densityContrast))
   
    # Generate the particle positions
    rotMatrix = randonRotationMatrix()
    r, p, pInSphere, t = particleLoop(nn0, rBoundary, halfWidthOfDomain, rotMatrix, Cmm, interpolatedRadius, ffB, fffB, tempFactor)
                                
    # Cull zero values
    pos = np.zeros((3, p))
    pos[0] = r[0][:p] * 3.09e18
    pos[1] = r[1][:p] * 3.09e18
    pos[2] = r[2][:p] * 3.09e18
    pTemp = t[:p]
                                
    # Work out mass of the particles
    ngas = p
    pMass = np.ones(p) * ((BEmass/2) * 1.991e33 / pInSphere) 

    return pos, pMass, ngas, pTemp, rBoundary

# Function to create a random rotation matrix to rotate particles by
def randonRotationMatrix():
    # Calculate sin and cos of psi
    psi = 2 * np.pi * np.random.random(1)
    cosPsi = np.cos(psi)
    sinPsi = np.sin(psi)
    
    # Calculate sin and cos of theta
    theta = np.random.random(1)
    cosTheta = 2 * theta - 1
    sinTheta = np.sqrt(1 - cosTheta**2)
    
    # Calculate sin and cos of phi
    phi = 2 * np.pi * np.random.random(1)
    cosPhi = np.cos(phi)
    sinPhi = np.sin(phi)
    
    # Create the matrix
    rotationMatrix = np.zeros((3,3))
    
    # Assign values of the matrix
    rotationMatrix[0,0] = (cosPhi * cosPsi) - (sinPhi * cosTheta * sinPsi)
    rotationMatrix[1,0] = (cosPhi * sinPsi) + (sinPhi * cosTheta * cosPsi)
    rotationMatrix[2,0] = sinPhi * sinTheta
    rotationMatrix[0,1] = - (sinPhi * cosPsi) - (cosPhi * cosTheta * sinPsi)
    rotationMatrix[1,1] = - (sinPhi * sinPsi) + (cosPhi * cosTheta * cosPsi)
    rotationMatrix[2,1] = cosPhi * sinTheta 
    rotationMatrix[0,2] = sinTheta * sinPsi
    rotationMatrix[1,2] = - sinTheta * cosPsi
    rotationMatrix[2,2] = cosTheta  
    
    return rotationMatrix 

# Function to adjust the properties of the BE sphere to make it collapse
def adjustProperties(pos, vels, pMass, rBoundary):
    # Calculate centre of mass and radial distance to it
    com = np.sum(pos[0] * pMass) / np.sum(pMass)
    rCentre = np.sqrt((pos[0] - com)**2 + (pos[1] - com)**2 + (pos[2] - com)**2)
    
    # Find padding particles
    padding = np.where(rCentre >= rBoundary*3.09e18)
    inside = np.where(rCentre < rBoundary*3.09e18)
    
    # Reset their velocites
    vels[0][padding] = 0
    vels[1][padding] = 0
    vels[2][padding] = 0
    
    # Adjust the masses
    pMass[padding] *= 0.1
    pMass[inside] *= 2
    
    return vels, pMass
    