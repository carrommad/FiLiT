# FiLiT
This function traces magnetic field lines to compute Poincaré plots and connection length in tokamaks and stellarators given its 3-dimensional magnetic field components. 

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


### 1.3. Euler vs. Runge-Kutta


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
  
### 4.1 filit.py
#### Configuration
Configure the run prior to execution:
  - select target mode between the following: 'matrix', 'line-horizontal', 'line-vertical'
    > tar_mode = 'line-horizontal'
  - choose number of targets per line:
    > ntargets1D = 100
  - choose number of toroidal turns:
    > nturns = 50
  - select the step size in mm of the field line for interpolation (this modifies variable *n_grid_phi_enhanced*; you can do it _ad hoc_):
    > step_phi_mm = 1
  - write the name (full path) of magnetic equilibrium 
    > input_file   = 'snapfile.80.nc'
  - full limiter limiter path (optional); if none, write the following:
    > lim_path     = ''
    
#### Execution
To run the file do:
> python3 filit.py
 

### 4.2. fiload.py


## Branches
### v1.2 (default)
- All interpolation functions are Python imported.
- First loop on turns and then loop on particles.
- Step size definition **not uniform** for every point but good approximation.
- Plotting improvements were carried out, still work to do.
- Difficult to parallelize but **fastest version**.
- Possible **change of magnetic field periodicity** with respect to the limiter.
- 2D geometries not implemented but easy.

### v1.1 (old)
- All interpolation functions are Python imported.
- First loop on turns and then loop on particles.
- Step size definition **not uniform** for every point but good approximation.
- Plotting improvements were carried out, still work to do.
- Difficult to parallelize but **fastest version**.
- Possibility to change magnetic field periodicity **not implemented**.

### v1.0 (old)
- All interpolation functions are manually written.
- First loop on turns and then loop on particles.
- Step size definition **not uniform** for every point but good approximation.
- Plotting improvements were strongly needed.
- Possibility to change magnetic field periodicity **not implemented**.

### v2.0 (slow)
- All interpolation functions are Python imported.
- First loop on particles and then loop on turns.
- Step size definition **uniform** for every point but good approximation.
- Plotting improvements were carried out, still work to do.
- Easier to parallelize particle loop but **slowest version**.
- Possibility to change magnetic field periodicity **not implemented**.

## Authors
- Carlos Romero
