#------------------------------------------------------------------------------
# MODULE - field_line_mod.py
#------------------------------------------------------------------------------
#> 
#> DESCRIPTION:
#> Contains auxiliary functions to trace magnetic field lines
#>
#> REVISION HISTORY:
#> 03 03 2021 - Carlos Romero
#>
#> TO DO:
#> -
#> 
#------------------------------------------------------------------------------
import numpy as np
#import IPython

#------------------------------------------------------------------------------
# FUNCTION - read_data_hint
#------------------------------------------------------------------------------
#> 
#> DESCRIPTION:
#> Reads netcdf output data from HINT snapfiles.
#> Adding or reducing returned variables is possible.
#> 
#> INPUT/OUTPUT:
#> - [IN]  dataset - netcdf data structure of an output (snapfile) from HINT
#> - [OUT] mtor    - toroidal mode number
#> - [OUT] rminb   - minimum value of R for the grid
#> - [OUT] rmaxb   - maximum value of R for the grid
#> - [OUT] zminb   - minimum value of Z for the grid
#> - [OUT] zmaxb   - maximum value of Z for the grid
#> - [OUT] B_R     - 3-tensor of radial magnetic field in the reactor
#> - [OUT] B_phi   - 3-tensor of toroidal magnetic field in the reactor
#> - [OUT] B_Z     - 3-tensor of vertical magnetic field in the reactor
#> 
#> EXAMPLE:
#> - Example:
#>      ds = nc.Dataset('snapfile.20.nc')
#>      rminb = read_data_hint(ds)[0]
#>      $ rminb = 2.52
#>
#> REVISION HISTORY:
#> 03 03 2021 - Carlos Romero
#> 
#> TO DO:
#> -
#> 
#------------------------------------------------------------------------------
def read_data_hint(dataset):

    mtor  = dataset.variables['mtor'][:]
    kstep = dataset.variables['kstep'][:]

    time  = dataset.variables['time'][:]
    rminb = dataset.variables['rminb'][:]
    rmaxb = dataset.variables['rmaxb'][:]
    zminb = dataset.variables['zminb'][:]
    zmaxb = dataset.variables['zmaxb'][:]

    B_R   = dataset.variables['B_R'][:]
    B_phi = dataset.variables['B_phi'][:]
    B_Z   = dataset.variables['B_Z'][:]

    v_R   = dataset.variables['v_R'][:]
    v_phi = dataset.variables['v_phi'][:]
    v_Z   = dataset.variables['v_Z'][:]

    P     = dataset.variables['P'][:]

    return mtor, rminb, rmaxb, zminb, zmaxb, B_R, B_phi, B_Z 

#------------------------------------------------------------------------------
# FUNCTION - create_targets
#------------------------------------------------------------------------------
#> 
#> DESCRIPTION:
#> Creates starting target points provided boundaries a number of points in 
#> one dimension. If 'matrix' mode selected, total targets equals ntargets^2.
#> 
#> INPUT/OUTPUT:
#> - [IN]  mode      - distribution of target particles: 'line' or 'matrix'
#> - [IN]  ntargets  - number of target particles of one dimension
#> - [IN]  rmin      - minimum radial coordinate
#> - [IN]  rmax      - maximum radial coordinate
#> - [IN]  zmin      - minimum vertical coordinate
#> - [IN]  zmax      - maximum vertical coordinate
#> - [OUT] targets_R - R-coordinates of target particles
#> - [OUT] targets_Z - Z-coordinates of target particles
#> 
#> EXAMPLE:
#> - Example:
#>
#> REVISION HISTORY:
#> 17 02 2021 - Carlos Romero
#> 
#> TO DO:
#> - Create example
#> -
#> 
#------------------------------------------------------------------------------
def create_targets(mode, ntargets, rmin, rmax, zmin, zmax):

    if mode == 'line':
        targets_R = np.linspace(rmin,rmax,ntargets)
        print('\n Default line choice: Z = (zmax - zmin) / 2 \n')
        targets_Z = (zmax - zmin) /2

    elif mode == 'matrix':
        targets_R = np.linspace(rmin,rmax,ntargets)
        targets_Z = np.linspace(zmin,zmax,ntargets)

    targets_R, targets_Z = np.meshgrid(targets_R, targets_Z)

    targets_R = targets_R.flatten()
    targets_Z = targets_Z.flatten()
    
    return targets_R, targets_Z

#------------------------------------------------------------------------------
# FUNCTION - grid_neighbours
#------------------------------------------------------------------------------
#> 
#> DESCRIPTION:
#> Given an array of coordinates (R,Z) and a grid of R-Z space, computes the
#> indices of grid for which the coordinates are closest and the coordinate 
#> value is just above the grid's value  
#> 
#> INPUT/OUTPUT:
#> - [IN]  targets_R - R-coordinates of target points
#> - [IN]  targets_Z - Z-coordinates of target points
#> - [IN]  grid_R    - R-coordinates of grid points
#> - [IN]  grid_Z    - Z-coordinates of grid points
#> - [OUT] index_R   - R-indices of grid points closest but below target point 
#> - [OUT] index_Z   - Z-indices of grid points closest but below target point 
#> - [OUT] dist_R    - distance from target point to closest but smallergrid
#>                   point in R
#> - [OUT] dist_Z    - distance from target point to closest but smallergrid
#>                   point in Z
#> 
#> EXAMPLE:
#> - Example:
#>      grid_R   = [0.  0.25 0.5  0.75 1.]
#>      grid_Z   = [0.  0.5  1.   1.5  2.]
#>      targets_R = [0.4 0.8]
#>      targets_Z = [0.  0.4]
#>      index_R, index_Z, dist_R, dist_Z = &
#>                          grid_neighbours(targets_R, targets_Z, grid_R, grid_Z)
#>      $ index_R = [1 3]
#>      $ index_Z = [1 0]
#>      $ dist_R  = [0.15 0.05]
#>      $ dist_Z  = [0.   0.4 ]
#>      $ grid_R[index_R] = [0.25 0.75]
#>      $ grid_Z[index_Z] = [0.   0.  ]
#> 
#> REVISION HISTORY:
#> 03 03 2021 - Carlos Romero
#> 
#> TO DO:
#> - 
#> 
#------------------------------------------------------------------------------
def grid_neighbours(targets_R, targets_Z, grid_R, grid_Z):

	index_R = np.zeros(targets_R.shape)
	index_Z = np.empty(targets_Z.shape)
	
	dist_R  = np.empty(targets_R.shape)
	dist_Z  = np.empty(targets_Z.shape)
	
	for m in range(len(targets_R)):
	
		r = np.argmin(abs(targets_R[m] - grid_R))
		if (targets_R[m] - grid_R[r] < 0):
			r = r - 1
			
		z = np.argmin(abs(targets_Z[m] - grid_Z))
		if (targets_Z[m] - grid_Z[z] < 0):
			z = z - 1
		
		index_R[m] = r
		index_Z[m] = z
		dist_R[m] = targets_R[m] - grid_R[r]
		dist_Z[m] = targets_Z[m] - grid_Z[z]

	return index_R, index_Z, dist_R, dist_Z

#------------------------------------------------------------------------------
# FUNCTION - interpolate_field
#------------------------------------------------------------------------------
#> 
#> DESCRIPTION:
#> Provided indices of a 2-dimensional grid in space, previously computed
#> distances of target points to points in that grid, and a tangent field,
#> computes for an interpolation of the field in the target point, which is not
#> necessarily in the grid.
#> 
#> INPUT/OUTPUT:
#> - [IN]  index_R - R-indices of grid points closest but below target
#> - [IN]  index_Z - Z-indices of grid points closest but below target
#> - [IN]  dist_R  - distance from target point to closest but smaller grid
#>                 point in R
#> - [IN]  dist_Z  - distance from target point to closest but smaller grid
#>                 point in Z
#> - [IN]  F_grid  - original scalar field matrix
#> - [IN]  step_R  - distance between two grid points in the R-coordinate
#> - [IN]  step_Z  - distance between two grid points in the Z-coordinate
#> - [OUT] F_inter - interpolated scalar field for targets
#> 
#> EXAMPLE:
#> - Example:
#>     index_R = np.array([0.])
#>     index_Z = np.array([0.])
#>     dist_R = np.array([0.])
#>     dist_Z = np.array([0.5])
#>     
#>     F_grid = np.empty((2,2))
#>     F_grid[0,:] = np.array([1.0, 2.0])
#>     F_grid[1,:] = np.array([ 3.0, 4.0])
#>     
#>     $ F_inter = [2.]
#>     
#>     
#> REVISION HISTORY:
#> 03 03 2021 - Carlos Romero
#> 
#> TO DO:
#> - 
#> 
#------------------------------------------------------------------------------
def interpolate_field(index_R, index_Z, dist_R, dist_Z, F_grid, step_R, step_Z):

	F_inter = np.empty(index_R.shape)
	
	for m in range(len(index_R)):
		r = int(index_R[m])
		z = int(index_Z[m])
		
		# Careful with the order of the coordinates [Z,R]
		F_inter[m] = (
					  F_grid[z,r]     * (step_R - dist_R[m]) * (step_Z - dist_Z[m]) + \
					  F_grid[z,r+1]   * (step_R - dist_R[m]) *           dist_Z[m]  + \
					  F_grid[z+1,r]   *           dist_R[m]  * (step_Z - dist_Z[m]) + \
					  F_grid[z+1,r+1] *           dist_R[m]  *           dist_Z[m]    \
					 )                                                                \
					 / (step_R * step_Z)
		
	return F_inter
	
#------------------------------------------------------------------------------
# FUNCTION - step_phi
#------------------------------------------------------------------------------
#> 
#> DESCRIPTION:
#> Traces discretized field lines in cylindrical geometry: field line equation.
#> CAREFUL with the conversion from toroidal to cyllindrical geometry!
#> 
#> INPUT/OUTPUT:
#> - [IN]    B_R_inter   - interpolated radial magnetic field for targets
#> - [IN]    B_phi_inter - interpolated toroidal magnetic field for targets
#> - [IN]    B_Z_inter   - interpolated vertical magnetic field for targets
#> - [IN]    delta_phi   - increment in toroidal direction
#> - [IN]    rmin        - minimum radial coordinate
#> - [IN]    rmax        - maximum radial coordinate
#> - [INOUT] targets_R   - R-coordinates of target points
#> - [INOUT] targets_Z   - Z-coordinates of target points
#> 
#> EXAMPLE:
#> - Example: unnecessary
#>
#> REVISION HISTORY:
#> 03 03 2021 - Carlos Romero
#> 
#> TO DO:
#> -
#> 
#------------------------------------------------------------------------------
def step_phi(B_R_inter, B_phi_inter, B_Z_inter, delta_phi, rmin, rmax, targets_R, targets_Z):
		
	# - the order is crucial: first update the Z coordinate!
	targets_Z = targets_Z + delta_phi * ( B_Z_inter / B_phi_inter ) * targets_R
	targets_R = targets_R + delta_phi * ( B_R_inter / B_phi_inter ) * targets_R
		
	return targets_R, targets_Z

#------------------------------------------------------------------------------
# FUNCTION - delete_outer_particles
#------------------------------------------------------------------------------
#> 
#> DESCRIPTION:
#> Deletes particles that are outside simulation boundaries.
#> 
#> INPUT/OUTPUT:
#> - [IN]    rmin      - minimum radial coordinate
#> - [IN]    rmax      - maximum radial coordinate
#> - [IN]    zmin      - minimum vertical coordinate
#> - [IN]    zmax      - maximum vertical coordinate
#> - [IN]    step_R    - distance between two grid points in the R-coordinate
#> - [IN]    step_Z    - distance between two grid points in the Z-coordinate
#> - [INOUT] targets_R - R-coordinates of target points
#> - [INOUT] targets_Z - Z-coordinates of target points
#> 
#> EXAMPLE:
#> - Example: omitting unnecessary variables for illustration
#>      rmin = 0
#>      targets_R = [2. 8.3 -1.3 2.7]
#>      targets_R = delete_outer_particles(targets_R, ~, rmin, ~, ~, ~)
#>      $ targets_R = [2. 8.3 2.7]
#>
#> REVISION HISTORY:
#> 24 02 2021 - Carlos Romero
#> 
#> TO DO:
#> - 
#> 
#------------------------------------------------------------------------------
def delete_outer_particles(rmin, rmax, zmin, zmax, step_R, step_Z, targets_R, targets_Z):

	condition_R = (rmin <= targets_R) & (targets_R <= rmax)	
	condition_Z = (zmin <= targets_Z) & (targets_Z <= zmax)
	
	targets_R   = targets_R[condition_R & condition_Z]
	targets_Z   = targets_Z[condition_Z & condition_R]
	
	return targets_R, targets_Z

#------------------------------------------------------------------------------
# FUNCTION - save_intersection
#------------------------------------------------------------------------------
#> 
#> DESCRIPTION:
#> Appends target vectors only when a turn has been completed (see conditional)
#> in loop. 
#> 
#> INPUT/OUTPUT:
#> - [IN]    targets_R   - R-coordinates of target points
#> - [IN]    targets_Z   - Z-coordinates of target points
#> - [INOUT] intersect_R - appended R-coordinates of target points at every turn
#> - [INOUT] intersect_Z - appended Z-coordinates of target points at every turn
#> 
#> EXAMPLE:
#> - Example: unnecessary
#>
#> REVISION HISTORY:
#> 17 02 2021 - Carlos Romero
#> 
#> TO DO:
#> - 
#> 
#------------------------------------------------------------------------------
def save_intersection(targets_R, targets_Z, intersect_R, intersect_Z):

    intersect_R = np.append(intersect_R,targets_R)
    intersect_Z = np.append(intersect_Z,targets_Z)

    return intersect_R, intersect_Z

#------------------------------------------------------------------------------
