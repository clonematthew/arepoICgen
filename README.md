[![PyPI package](https://img.shields.io/badge/pip%20install-arepoICgen-brightgreen)](https://pypi.org/project/arepoICgen/) 
[![version number](https://img.shields.io/pypi/v/arepoICgen?color=green&label=version)](https://github.com/clonematthew/arepoICgen/releases)

## arepoICgen

Package for creating Initial Conditions for the AMR code AREPO[^1] for simulations of molecular clouds intended for scientific research. Many parts of this code are adapted from a FORTRAN code developed by Paul Clark called _setpartarepo_. 

This code is available via a direct clone of this repository, or through pip: `pip install arepoICgen`

#### Usage Guide

arepoICgen takes in two python dictionaries as inputs, **config** and **params**. Config provides the code information about the geometry and metadata of the cloud, while params supplies the physical parameters of the cloud. Some config and param flags are required, and others are optional. Required flags will be marked with a "*".

##### Config

| Flag | Options | Description | Associated Params |
| ----------- | ----------- | ---- | --| 
| verbose | True, False | Whether to print extended feedback as the code runs (does not affect the ICs) | |
| grid* | boxGrid, boxRan, sphereGrid, sphereRan, ellipseRan, cylinderRan | The geometry and particle distribution (gridded or random). | lengths, radii |
| turbulence | turbFile  | Whether to apply to turbulent velocity field from a file | virialParam | 
| turbFile | _path_to_file_ | The path to a turbulent velocity grid file, needed only for _turbFile_ | | 
| turbSize | _integer_ | The size of the velocity grid (assumed equal length cubic), needed only for _turbFile_ | |
| rotation | rotation| Whether to add rotation to the cloud | beta, rotationRadius |
| extras | bossBodenheimer, densityGradient, bonnorEbert | What density perturbation to add to the cloud, either Boss Bodenheimer[^2], a linear density gradient[^3] or a Bonnor-Ebert[^4] density profile (note this setting automatically outputs density as mass) |
| padding* | True, False | Whether to pad the cloud with low density particles | tempFactor, boxDims, paddingDensity |
| bField | True, False | Whether we need to give the particles an inital magnetic field value (note this does not do anything besides gives particles an initial B field of 0) | |
| outValue | masses, density | Option to output the mass field as density (note this needs AREPO to be compiled with the right config flag for this) | density |
| output* | hdf5 | Filetype to write the ICs as (type 2 binary coming... one day) | |
| filename* | "..." | Name of the IC file (note that the file extension is not needed) | | 

##### Params

Note that lengths and radii are both marked as required, but only one may be needed in some cases. 

| Flag | Format/Units | Description | Associated Config |
| ---- | ------ | ----------- | ----------------- |
| ngas* |  | Number of particles to populate the cloud with | |
| lengths* | [xLength, yLength, zLength] $\rm [pc]$ | The length of the box or ellipse in each direction (for the cylinder, xLength is taken as the length of the cylinder) | grid |
| radii* | Radius $\rm [pc]$ | The radius of a spherical cloud, or the radius of the face of a cylinder (not needed by box or elliptical setups) | grid | 
| mass* | $\rm [M_\odot]$ | The total mass of the cloud | |
| temp* | $\rm [K]$ | The inital temperature of the gas inside the cloud | |
| mu* | $\mu$ | The mean molecular weight of the gas in the cloud, typically 1.4 for atomic and 2.4 for molecular | |
| virialParam | $\alpha = \rm\frac{E_{Kin}}{E_{Grav}}$ | The virial parameter of the cloud, used for scaling turbulent velocities | turbulence |
| beta | $\beta = \rm \frac{E_{Rot}}{E_{Grav}}$ | The ratio between rotational and gravitational energy | rotation |
| rotationRadius | AU | The radius outside of which to add rotation. Allows no angular momentum in the centre of the cloud so that numerical rings do not form | rotation |
| boxSize| [x, y, z] | The factor by which the box should be larger than the cloud in each dimension (i.e 2x, 3x, 4x the size) | padding |
| tempFactor | | How much hotter to make the padding particles compared to the cloud particles | padding |
| paddingDensity |   | What fraction of the cloud's density should the padding region be (ignored for Bonnor-Ebert profile, padding density set to match outer BE sphere density) | padding | 
| density | $\rm cm^{-3}$ | Density to give particles when outputting mass as density | outValue |

#### Calling generateICs

arepoICgen is called using code similar to the following, where it is wrapped in the `generateICs` function:

```
from arepoICgen import generateICs

config = {}
params = {}

generateICs(config, params)
```

##### Examples

See examples.txt for some example setups of the config and params dictionaries. 

[^1]: AREPO is publicly available, see [here](https://arepo-code.org)
[^2]: See Boss & Bodenheimer (1979)
[^3]: See Bonnell et al. (2008) 
[^4]: See Bonnor (1956) and Ebert (1955)
