# FieldLineTracer
This function traces magnetic field lines to compute Poincaré plots and connection length in tokamaks and stellarators given its 3-dimensional magnetic field components. 


# 1. Definition of Field Line
A magnetic field line is a 3D curve for which the magnetic field is tangent at every point.
Define $\vec{\Gamma} \in \mathbf{R}^3$ as our 3D curve parametrized by $s \in \mathbf{R}$:
$$ \frac{\partial \vec{\Gamma}}{\partial s} = \vec{B}$$

For simplication, consider the different component of our cuve as follows: $ \Gamma_R = R $, $ \Gamma_{\phi} = \phi $ and $ \Gamma_Z = Z $.

To obtain the corresponding scalar components of the differential equation, the Jacobian needs to be taken into account.


## 1.1. Cylindrical Coordinates
The differential equation in cylindrical coordinates $[R, \phi, Z]$ for the field lines being:

- $$ \frac{\partial R}{\partial \phi} =  R \cdot \frac{B_R}{B_{\phi}}$$
- $$ \frac{\partial Z}{\partial \phi} =  R \cdot \frac{B_Z}{B_{\phi}}$$


being the finite differences:

- $$ R(k+1) = R(k) + \Delta \phi \cdot  R(k) \cdot \frac{B_R(k)}{B_{\phi}(k)}$$
- $$ Z(k+1) = Z(k) + \Delta \phi \cdot  R(k) \cdot \frac{B_Z(k)}{B_{\phi}(k)}$$


## 1.2. Toroidal Coordinates
The differential equation in toroidal coordinates [R, \phi, \theta] for the field lines being:
**Pending implementation**


# 2. Limiter Geometry
If a limiter or divertor 3-dimensional geometry is available, an option can be selected to compute more realistic trajectories.


# 3. INPUTS
To run this program, the following files are necessary:
  - a 3D magnetic equilibrium netcdf file containing radial, toroidal and vertical components of the field
  - a 3D limiter netcdf file containing information on the inside/outside for the simulation
 
# 4. EXAMPLE
Links to example files:
  - magnetic equilibrium: https://mega.nz/file/EwsVUIiA#ZyX3waZSOhuUkpyO_aMJv_z_cmZKIGY9fqAY1HmpriM
  - limiter file: https://mega.nz/file/Y1dl1aCQ#_-FDM5VHMNyA5Q2acpxoYyqjTxWYKoSuEYIQ5iu2Wuc
  
# 5. RUNNING THE EXAMPLE
To run the example:
  - write the name of the files (full path): magnetic equilibrium and limiter
  - select if you want the limiter to be included with 'lim_mode'
  - select if you want to compute the connection length with 'lc_mode'
  - select if you want to print and plot debug options with 'debug'


# 6. Currently working version
At the moment, 'changes/Python' branch is the most updated. A pull request to 'main' branch of the *Python* version and completion of the *Fortran* version will be done soon.

## 6.1 2D version
2D magnetic field and limiter can be easily replicated in the toroidal dimension to create a 3D version.


# Authors
- Carlos Romero
