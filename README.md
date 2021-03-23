# FiLiT
This function traces field lines for any 3D vector field to compute Poincaré plots and connection length. Particularly, for tokamaks and stellarators, given its 3-dimensional magnetic field components. 

## 1. Definition of Field Line
A field line for a given vector field is a 3D curve for which the vector field is tangent at every point in space.

We use a cylindrical reference frame for the vector fields *(R, phi, Z)* to follow field lines but, thanks to the modularity of the program, this can be easily modified to other geometries

## 2. Limiter Geometry
If a limiter geometry is available, the option can be selected to compute more realistic trajectories, dynamically deleting particles that go outside our limiter boundaries.

## 3. Inputs
To run this program, the following files are necessary:
  - **3D cylindrical (magnetic) field** (for magnetic fields, preferably an equilibrium) *netcdf* file containing radial, toroidal and vertical components of the field:
    > https://mega.nz/file/EwsVUIiA#ZyX3waZSOhuUkpyO_aMJv_z_cmZKIGY9fqAY1HmpriM
  - **3D limiter (optional)** *netcdf* file containing information on the inside/outside for the simulation:
    > https://mega.nz/file/Y1dl1aCQ#_-FDM5VHMNyA5Q2acpxoYyqjTxWYKoSuEYIQ5iu2Wuc
    
## 4. Branches
### v1.2 (default)
- All interpolation functions are Python imported.
- First loop on turns and then loop on particles.
- Step size definition **not uniform** for every point but good approximation.
- Plotting improvements were carried out, still work to do.
- Difficult to parallelize but **fastest version**.
- Possible **change of magnetic field periodicity** with respect to the limiter.
- 2D geometries not implemented but easy.
- Contains **Fortran** version

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
- Carlos Romero Madrid
