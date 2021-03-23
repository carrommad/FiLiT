## 1. filit.F90
This is the main field line tracing **Fortran** program.
### 1.1. Configuration
Configure the run prior to execution:
  - select target mode between the following: 'matrix', 'line-horizontal', 'line-vertical'
    > tar_mode = 'line-horizontal'
  - choose number of targets per line:
    > ntargets1D = 100
  - choose number of toroidal turns:
    > nturns = 50
  - select the step size in mm of the field line for interpolation (this modifies variable **n_grid_phi_enhanced**; you can do it *ad hoc*):
    > step_phi_mm = 1
  - write the name (full path) of magnetic equilibrium 
    > input_file   = 'snapfile.80.nc'
  - full limiter limiter path (optional); if none, write the following:
    > lim_path     = ''
    
### 1.2. Execution
To run the script do:
  > (***???***)
 
Depending on the number of targets, number of turns and plotting resolution (see *filib.f90*), it can take from 1-100 (***???***) minutes to compute and plot a field line tracing experiment.
 

## 2. fiload.f90
This is an auxiliary **Fortran** script to load and plot the results of the main program.
### 2.1. Configuration
- full path of *filit.f90* output file:
  > filename = 'results/matrix_ntargets-10000_nturns-500'
- boolean for plotting initial target distribution:
  > plot_initial = 0
  
### 2.2. Execution
To run the script do:
  > (***???***)

Depending on the number of targets, number of turns and plotting resolution (see *filib.f90*), it can take from 1-10 (***???***) minutes to load and plot a *filit.f90* output file.

## 3. filib.py
This is an auxiliary **Fortran** script containing all the functions to run *filit.f90* and *fiload.f90*:
-  **read_data_hint**: reads magnetic equilibrium data from *HINT* output file (see https://github.com/yasuhiro-suzuki/HINT3D.git).
-  **read_limiter**: reads limiter data from *MKLIM* output file (see https://github.com/yasuhiro-suzuki/HINT3D.git).
-  **limiter_boundary**: computes the boundaries of the limiter just for plotting.
-  **create_targets**: creates a target distribution in a line (horizontal or vertical) or a matrix. This function can be easily improved.
-  **step_phi_cylindrical_rk4**: 4th order Runge-Kutta specific for the following field line cylindrical equation:
    - **field_line_function**: cylindrical equation of field line to be evaluated at the differents steps of a Runge-Kutta scheme.
-  **delete_outer_particles**: deletes particles outside computational boundaries.
-  **delete_outer_particles_limiter**: deletes particles outside dynamical limiter.
-  **plot_Poincare_Lc**: this representation script provides Poincaré plots and connection length plots for the field line tracing done to our magnetic equilibrium.
