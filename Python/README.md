## Example
Links to example files:
  - 3D magnetic equilibrium in netcdf format: https://mega.nz/file/EwsVUIiA#ZyX3waZSOhuUkpyO_aMJv_z_cmZKIGY9fqAY1HmpriM
  - 3D limiter file in netcdf format: https://mega.nz/file/Y1dl1aCQ#_-FDM5VHMNyA5Q2acpxoYyqjTxWYKoSuEYIQ5iu2Wuc
  
## filit.py
This is the main field line tracing Python program.
### Configuration
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
    
### Execution
To run the script do:
> python3 filit.py
 

## fiload.py
This is an auxiliary Python script to load and plot the results of the main program.
### Configuration
- full path of *filit.py* output file:
  > filename = 'results/matrix_ntargets-10000_nturns-500'
- boolean for plotting initial target distribution:
  > plot_initial = 0
  
### Execution
To run the script do:
> python3 fiload.py


## filib.py
This is an auxiliary Python script containing all the functions to run *filit.py* and *fiload.py*:
-  **read_data_hint**: reads magnetic equilibrium data from *HINT* output file (see https://github.com/yasuhiro-suzuki/HINT3D.git).
-  **read_limiter**: reads limiter data from *MKLIM* output file (see https://github.com/yasuhiro-suzuki/HINT3D.git).
-  **limiter_boundary**: computes the boundaries of the limiter just for plotting.
-  **create_targets**: creates a target distribution in a line (horizontal or vertical) or a matrix. This function can be easily improved.
-  **step_phi_cylindrical_rk4**: 4th order Runge-Kutta specific for the following field line cylindrical equation:
  - ** field_line_function**: cylindrical equation of field line to be evaluated at the differents steps of a Runge-Kutta scheme.
-  **delete_outer_particles**: deletes particles outside computational boundaries.
-  **delete_outer_particles_limiter**: deletes particles outside dynamical limiter.
-  **plot_Poincare_Lc**: this representation script provides Poincaré plots and connection length plots for the field line tracing done to our magnetic equilibrium.

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
