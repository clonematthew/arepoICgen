# Importing libraries
import numpy as np

# Function to add solid body rotation
def addRotation(pos, pMass, vels, beta):
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
    omega = np.sqrt( 6.673-8 * 3. * beta * mtot / (rMax**3))

    # Adding the rotation to the x and y velocities, rotating about the z axis
    vels[0] -= omega * (pos[1] - ycom)
    vels[1] += omega * (pos[0] - xcom)

    # Working out gravitational potential energy
    eGrav = (6.67e-8) * (3./5.) * (mtot**2) / rMax
    
    # Working out rotational kinetic energy
    r = (pos[0] - xcom)**2 + (pos[1] - ycom)**2
    vRot = (pos[0] - xcom) * vels[1] - (pos[1] - ycom) * vels[0]
    eRot = np.sum(0.5 * pMass * vRot * vRot / r)

    # Reporting the deviation from the desired beta value
    print("Difference from desired beta: {:.1f}%".format(abs(100*(beta-eRot/eGrav)/beta)))

    return vels