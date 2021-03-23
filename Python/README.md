## Example
Links to example files:
  - 3D magnetic equilibrium in netcdf format: https://mega.nz/file/EwsVUIiA#ZyX3waZSOhuUkpyO_aMJv_z_cmZKIGY9fqAY1HmpriM
  - 3D limiter file in netcdf format: https://mega.nz/file/Y1dl1aCQ#_-FDM5VHMNyA5Q2acpxoYyqjTxWYKoSuEYIQ5iu2Wuc
  
## filit.py
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
To run the file do:
> python3 filit.py
 

## fiload.py


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
