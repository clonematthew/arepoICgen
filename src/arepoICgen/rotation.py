# Importing libraries
import numpy as np

# Function to add solid body rotation
def addRotation(pos, pMass, vels, beta, rotationRadius, verbose=False):
    # Find the centre of mass of the body
    mtot = np.sum(pMass) 
    xcom = np.sum(pMass * pos[0]) / mtot
    ycom = np.sum(pMass * pos[1]) / mtot
    zcom = np.sum(pMass * pos[2]) / mtot

    # Working out rmax of the sphere
    xmax = np.max(pos[0]) - xcom
    ymax = np.max(pos[1]) - ycom
    zmax = np.max(pos[2]) - zcom
    rMax = np.max([xmax, ymax, zmax]) 

    # Working out rotational velocity
    omega = np.sqrt(6.673e-8 * 3. * beta * mtot / (rMax**3))
    
    # Working out each cells distance from the centre
    rDist = np.sqrt((pos[0]-xcom)**2 + (pos[1]-ycom)**2 + (pos[2]-zcom)**2)
    outsideRadius = np.where(rDist > rotationRadius*1.5e13)

    # Adding the rotation to the x and y velocities, rotating about the z axis
    vels[0][outsideRadius] -= omega * (pos[1][outsideRadius] - ycom)
    vels[1][outsideRadius] += omega * (pos[0][outsideRadius] - xcom)

    # Working out gravitational potential energy
    eGrav = (6.67e-8) * (3./5.) * (mtot**2) / rMax
 
    # Working out the rotational energy
    momentOfInteria = (2/5) * mtot * rMax**2
    eRot = (1/2) * momentOfInteria * omega**2

    if verbose:
        # Reporting the deviation from the desired beta value
        print("Difference from desired beta: {:.2f}%".format(abs(100*(beta-eRot/eGrav)/beta)))

    return vels