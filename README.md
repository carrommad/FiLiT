# FiLiT
This function traces magnetic field lines to compute Poincaré plots and connection length in tokamaks and stellarators given its 3-dimensional magnetic field components. 

## 1. Definition of Field Line
A magnetic field line is a 3D curve for which the magnetic field is tangent at every point.

We use a cylindrical magnetic field (R, phi, Z) to follow field lines but, thanks to the modularity of the program, this can be easily modified to other geometries

## 2. Limiter Geometry
If a limiter geometry is available, the option can be selected to compute more realistic trajectories, dynamically deleting particles that go outside our limiter boundaries.

## 3. Inputs
To run this program, the following files are necessary:
  - **3D cylindrical magnetic field** (preferably, an equilibrium) *netcdf* file containing radial, toroidal and vertical components of the field:
    > https://mega.nz/file/EwsVUIiA#ZyX3waZSOhuUkpyO_aMJv_z_cmZKIGY9fqAY1HmpriM
  - **3D limiter (optional)** *netcdf* file containing information on the inside/outside for the simulation:
    > https://mega.nz/file/Y1dl1aCQ#_-FDM5VHMNyA5Q2acpxoYyqjTxWYKoSuEYIQ5iu2Wuc

## Authors
- Carlos Romero Madrid
