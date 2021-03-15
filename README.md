# FieldLineTracer
This function traces magnetic field lines to compute Poincaré plots and connection length in tokamaks and stellarators given its 3-dimensional magnetic field components. 

## 0. Branches
- *v1.0: This is the primary version where all functions were manually written*
- ***v1.1: ***
- *v2.0:*

## 1. Definition of Field Line
A magnetic field line is a 3D curve for which the magnetic field is tangent at every point.
Define $\vec{\Gamma} \in \mathbf{R}^3$ as our 3D curve parametrized by $s \in \mathbf{R}$:
- $$ \frac{\partial \vec{\Gamma}}{\partial s} = \vec{B}$$

For simplication, consider the different component of our cuve as follows: 
- $ \Gamma_R = R $
- $ \Gamma_{\phi} = \phi $
- $ \Gamma_Z = Z $.

To obtain the corresponding scalar components of the differential equation, the Jacobian needs to be taken into account. 
- **Explain better and add Jacobian to each coordinate type**


### 1.1. Cylindrical Coordinates
The differential equation in cylindrical coordinates $[R, \phi, Z]$ for the field lines being:

- $$ \frac{\partial R}{\partial \phi} =  R \cdot \frac{B_R}{B_{\phi}}$$
- $$ \frac{\partial Z}{\partial \phi} =  R \cdot \frac{B_Z}{B_{\phi}}$$


considerin *k* as the toroidal index, the finite differences equation is:

- $$ R(k+1) = R(k) + \Delta \phi \cdot  R(k) \cdot \frac{B_R(k)}{B_{\phi}(k)}$$
- $$ Z(k+1) = Z(k) + \Delta \phi \cdot  R(k) \cdot \frac{B_Z(k)}{B_{\phi}(k)}$$


### 1.2. Toroidal Coordinates
The differential equation in toroidal coordinates $[R, \phi, \theta]$ for the field lines being:
- **Pending implementation**


## 2. Limiter Geometry
If a limiter or divertor 3-dimensional geometry is available, the option *lim_mod = 1* can be selected to compute more realistic trajectories.


## 3. Inputs
To run this program, the following files are necessary:
  - a 3D magnetic equilibrium netcdf file containing radial, toroidal and vertical components of the field
  - a 3D limiter netcdf file containing information on the inside/outside for the simulation
 
## 4. Example
Links to example files:
  - 3D magnetic equilibrium in netcdf format: https://mega.nz/file/EwsVUIiA#ZyX3waZSOhuUkpyO_aMJv_z_cmZKIGY9fqAY1HmpriM
  - 3D limiter file in netcdf format: https://mega.nz/file/Y1dl1aCQ#_-FDM5VHMNyA5Q2acpxoYyqjTxWYKoSuEYIQ5iu2Wuc
  
### 4.1 Running the Example
To run the example:
  - write the name of the files (full path): magnetic equilibrium and limiter
  - select if you want the limiter to be included with *lim_mode*
  - select if you want to compute the connection length with *lc_mode*
  - select if you want to print and plot debug options with *debug*
  - 
### 4.1 2D version
2D magnetic field and limiter can be easily replicated in the toroidal dimension to create a 3D version.

## Authors
- Carlos Romero
