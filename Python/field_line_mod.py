#==============================================================================
# MODULE - field_line_mod.py
#==============================================================================
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
#==============================================================================
import numpy as np
#import IPython

#==============================================================================
# FUNCTION - read_data_hint
#==============================================================================
#> 
#> DESCRIPTION:
#> Reads netcdf output data from HINT snapfiles.
#> Adding or reducing returned variables is possible.
#> 
#> INPUT/OUTPUT:
#> - [IN]  dataset - netcdf data structure of an output (snapfile) from HINT
#> - [OUT] mtor    - toroidal periodicity of reactor
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
#==============================================================================
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
    
#==============================================================================
# FUNCTION - read_limiter
#==============================================================================
#> 
#> DESCRIPTION:
#> Reads netcdf limiter/divertor geometry data.
#> 
#> INPUT/OUTPUT:
#> - [IN]  dataset - netcdf data structure of limiter/divertor geometry
#> - [OUT] limiter - 3-tensor with limiter/divertor geometry information
#> 
#> EXAMPLE:
#> - Example:
#>      ds = nc.Dataset('limiter.nc')
#>      limiter = read_limiter(ds)
#>      $ limiter[0,0,:] = [0, 0, ... 1, 1, 1, ... 0, 0]
#>
#> REVISION HISTORY:
#> 03 03 2021 - Carlos Romero
#> 
#> TO DO:
#> -
#> 
#==============================================================================
def read_limiter(dataset):
	
	limiter = dataset.variables['limiter'][:]
	
	return limiter
	
#==============================================================================
# FUNCTION - limiter_boundary
#==============================================================================
#> 
#> DESCRIPTION:
#> Given a toroidal slice with limiter/divertor information, computes the 
#> limiter boundary
#> 
#> INPUT/OUTPUT:
#> - [IN]  grid_R        - R-coordinates of grid points
#> - [IN]  grid_Z        - Z-coordinates of grid points
#> - [IN]  limiter_slice - limiter toroidal slice 
#> - [OUT] limiter_R     - R-coordinates of limiter line points
#> - [OUT] limiter_Z     - Z-coordinates of limiter line points
#> 
#> EXAMPLE:
#> - Example:
#>      ds = nc.Dataset('limiter.nc')
#>      limiter = read_limiter(ds)
#>      $ limiter[0,0,:] = [0, 0, ... 1, 1, 1, ... 0, 0]
#>
#> REVISION HISTORY:
#> 03 03 2021 - Carlos Romero
#> 
#> TO DO:
#> -
#> 
#==============================================================================
def limiter_boundary(grid_R, grid_Z, limiter_slice):
	
	# - compute slice dimensions
	Nz = limiter_slice.shape[0]
	Nr = limiter_slice.shape[1]
	
	# - declare lists to be appended
	limiter_R = []
	limiter_Z = []
	
	# - looping limiter values 
	for i in range(Nz):
		for j in range(Nr-1):
					
			# - stepping inside the limiter
			if limiter_slice[i,j+1] - limiter_slice[i,j] == 1:
			
				limiter_R.append(grid_R[j+1])
				limiter_Z.append(grid_Z[i])
				
			# - stepping outside the limiter
			elif limiter_slice[i,j+1] - limiter_slice[i,j] == -1:

				limiter_R.append(grid_R[j])
				limiter_Z.append(grid_Z[i])
			
	# - creating numpy arrays	
	limiter_R = np.array(limiter_R)
	limiter_Z = np.array(limiter_Z)
	
	return limiter_R, limiter_Z

#==============================================================================
# FUNCTION - create_targets
#==============================================================================
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
#==============================================================================
def create_targets(mode, ntargets, rmin, rmax, zmin, zmax):

    if mode == 'line':
        targets_R = np.linspace(rmin,rmax,ntargets)
        targets_Z = (zmax + zmin) / 2
        
        print('\n')
        print('=====================================================')
        print('  DEFAULT LINE CHOICE: Z = (zmax + zmin) / 2 = ' + str(targets_Z))
        print('=====================================================')
        print('\n')
        
    elif mode == 'matrix':
        targets_R = np.linspace(rmin,rmax,ntargets)
        targets_Z = np.linspace(zmin,zmax,ntargets)

    targets_R, targets_Z = np.meshgrid(targets_R, targets_Z)

    targets_R = targets_R.flatten()
    targets_Z = targets_Z.flatten()
    
    return targets_R, targets_Z

#==============================================================================
# FUNCTION - grid_neighbours_2D
#==============================================================================
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
#==============================================================================
def grid_neighbours_2D(targets_R, targets_Z, grid_R, grid_Z):

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
	
#==============================================================================
# FUNCTION - grid_neighbours_3D
#==============================================================================
#> 
#> DESCRIPTION:
#> Given an array of coordinates (R,Z), a toroidal angle (phi) and a tensor 
#> grid in [phi, Z, R] space, computes the indices of grid for which the
#> coordinates are closest and the coordinate value is just above the grid's
#> value  
#> 
#> INPUT/OUTPUT:
#> - [IN]  angle      - toroidal angle in radians of current iteration
#> - [IN]  n_grid_phi - number of slices in the original toroidal grid
#> - [IN]  targets_R  - R-coordinates of target points
#> - [IN]  targets_Z  - Z-coordinates of target points
#> - [IN]  grid_R     - R-coordinates of grid points
#> - [IN]  grid_Z     - Z-coordinates of grid points
#> - [OUT] index_R    - R-indices of grid points closest but below target point 
#> - [OUT] index_Z    - Z-indices of grid points closest but below target point 
#> - [OUT] index_phi  - phi-index of grid points closest but below target point 
#> - [OUT] dist_R     - distance from target point to closest but smallergrid
#>                    point in R
#> - [OUT] dist_Z     - distance from target point to closest but smallergrid
#>                    point in Z
#> - [OUT] dist_phi   - distance from target point to closest but smallergrid
#>                    point in phi
#> 
#> EXAMPLE:
#> - Example:
#>      grid_R     = [0.  0.25 0.5  0.75 1.]
#>      grid_Z     = [0.  0.5  1.   1.5  2.]
#>      targets_R  = [0.4 0.8]
#>      targets_Z  = [0.  0.4]
#>      angle      = 0.1
#>      n_grid_phi = 10
#>      index_R, index_Z, index_phi, dist_R, dist_Z, dist_phi = \
#>                          grid_neighbours(targets_R, targets_Z, grid_R, grid_Z)
#>      $ index_R   = [1 3]
#>      $ index_Z   = [1 0]
#>      $ index_phi = 0
#>      $ dist_R    = [0.15 0.05]
#>      $ dist_Z    = [0.   0.4 ]
#>      $ dist_phi  = 0.1
#>      $ grid_R[index_R]     = [0.25 0.75]
#>      $ grid_Z[index_Z]     = [0.   0.  ]
#>      $ grid_phi[index_phi] = 0.
#> 
#> REVISION HISTORY:
#> 09 03 2021 - Carlos Romero
#> 
#> TO DO:
#> - 
#> 
#==============================================================================
def grid_neighbours_3D(angle, n_grid_phi, targets_R, targets_Z, grid_R, grid_Z):

	index_R = np.zeros(targets_R.shape)
	index_Z = np.empty(targets_Z.shape)
	
	index_phi = int((n_grid_phi-1) * angle / ( 2 * np.pi ))
	
	dist_R  = np.empty(targets_R.shape)
	dist_Z  = np.empty(targets_Z.shape)
	
	dist_phi = angle - 2 * np.pi * index_phi / (n_grid_phi - 1)
	
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

	return index_R, index_Z, index_phi, dist_R, dist_Z, dist_phi

#==============================================================================
# FUNCTION - interpolate_field_2D
#==============================================================================
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
#==============================================================================
def interpolate_field_2D(index_R, index_Z, dist_R, dist_Z, F_grid, step_R, step_Z):

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
	
#==============================================================================
# FUNCTION - interpolate_field_3D
#==============================================================================
#> 
#> DESCRIPTION:
#> Provided indices of a 2-dimensional grid in space, previously computed
#> distances of target points to points in that grid, and a tangent field,
#> computes for an interpolation of the field in the target point, which is not
#> necessarily in the grid.
#> 
#> INPUT/OUTPUT:
#> - [IN]  index_R  - R-indices of grid points closest but below target
#> - [IN]  index_Z  - Z-indices of grid points closest but below target
#> - [IN]  index_phi  - phi-indices of grid points closest but below target
#> - [IN]  dist_R   - distance from target point to closest but smaller grid
#>                  point in R
#> - [IN]  dist_Z   - distance from target point to closest but smaller grid
#>                  point in Z
#> - [IN]  dist_phi   - distance from target point to closest but smaller grid
#>                  point in phi
#> - [IN]  F_grid   - original scalar field 3-tensor
#> - [IN]  step_R   - distance between two grid points in the R-coordinate
#> - [IN]  step_Z   - distance between two grid points in the Z-coordinate
#> - [IN]  step_phi - distance between two grid points in the phi-coordinate
#> - [OUT] F_inter  - interpolated scalar field for targets
#> 
#> EXAMPLE:
#> - Example:

#>     
#>     
#> REVISION HISTORY:
#> 09 03 2021 - Carlos Romero
#> 
#> TO DO:
#> - 
#> 
#==============================================================================
def interpolate_field_3D(index_R, index_Z, index_phi, dist_R, dist_Z, dist_phi, F_grid, step_R, step_Z, step_phi):

	F_inter = np.empty(index_R.shape)
	
	if index_phi == 63:
		index_next = 0
	else:
		index_next = index_phi + 1
	
	for m in range(len(index_R)):
		r = int(index_R[m])
		z = int(index_Z[m])
		
		# Careful with the order of the coordinates [Z,R]
		F_inter[m] = (                                                                               \
					   (                                                                             \
					   	 F_grid[index_phi,z,r]     * (step_R - dist_R[m]) * (step_Z - dist_Z[m]) +   \
					     F_grid[index_phi,z,r+1]   * (step_R - dist_R[m]) *           dist_Z[m]  +   \
					     F_grid[index_phi,z+1,r]   *           dist_R[m]  * (step_Z - dist_Z[m]) +   \
					     F_grid[index_phi,z+1,r+1] *           dist_R[m]  *           dist_Z[m]      \
					   )                                                                             \
					   * (step_phi - dist_phi)                                                       \
					 +                                                                               \
					   (                                                                             \
					   	 F_grid[index_next,z,r]     * (step_R - dist_R[m]) * (step_Z - dist_Z[m]) + \
					     F_grid[index_next,z,r+1]   * (step_R - dist_R[m]) *           dist_Z[m]  + \
					     F_grid[index_next,z+1,r]   *           dist_R[m]  * (step_Z - dist_Z[m]) + \
					     F_grid[index_next,z+1,r+1] *           dist_R[m]  *           dist_Z[m]    \
					   )                                                                             \
					   *                  dist_phi                                                   \
					 )                                                                               \
					 / (step_R * step_Z * step_phi)
		
	return F_inter
	
#==============================================================================
# FUNCTION - step_phi_cylindrical
#==============================================================================
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
#==============================================================================
def step_phi_cylindrical(B_R_inter, B_phi_inter, B_Z_inter, delta_phi, rmin, rmax, targets_R, targets_Z):
		
	# - the order is crucial: first update the Z coordinate!
	targets_Z = targets_Z + delta_phi * ( B_Z_inter / B_phi_inter ) * targets_R
	targets_R = targets_R + delta_phi * ( B_R_inter / B_phi_inter ) * targets_R
		
	return targets_R, targets_Z
	
#==============================================================================
# FUNCTION - delete_outer_particles
#==============================================================================
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
#> - [INOUT] active    - array containing information (0/1) regarding if 
#>                     particle is still inside boundaries
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
#==============================================================================
def delete_outer_particles(rmin, rmax, zmin, zmax, step_R, step_Z, targets_R, targets_Z, active):

	condition_R = (rmin <= targets_R) & (targets_R <= rmax)	
	condition_Z = (zmin <= targets_Z) & (targets_Z <= zmax)
	
	condition   = condition_R & condition_Z
	
	# - changing values of active particles only
	active[active == 1.] = condition
	
	targets_R   = targets_R[condition]
	targets_Z   = targets_Z[condition]
	
	return targets_R, targets_Z, active

#==============================================================================
# FUNCTION - delete_outer_particles_limiter
#==============================================================================
#> 
#> DESCRIPTION:
#> Deletes particles that are outside simulation boundaries.
#> 
#> INPUT/OUTPUT:
#> - [IN]    limiter_inter - interpolated limiter information of targets
#> - [INOUT] targets_R     - R-coordinates of target points
#> - [INOUT] targets_Z     - Z-coordinates of target points
#> - [INOUT] active        - array containing information (0/1) regarding if 
#>                         particle is still inside limiter/divertor geometry
#> 
#> EXAMPLE:
#> - Example: omitting unnecessary variables for illustration
#>      rmin = 0
#>      targets_R = [2. 8.3 -1.3 2.7]
#>      targets_R = delete_outer_particles(targets_R, ~, rmin, ~, ~, ~)
#>      $ targets_R = [2. 8.3 2.7]
#>
#> REVISION HISTORY:
#> 03 03 2021 - Carlos Romero
#> 
#> TO DO:
#> - 
#> 
#==============================================================================
def delete_outer_particles_limiter(limiter_inter, targets_R, targets_Z, active):

	# - if the particle is close to the wall, its interpolation value is
	# between 0 and 1
	condition = limiter_inter != 0.
	
	# - changing values of active particles only
	active[active == 1.] = condition
	
	targets_R   = targets_R[condition]
	targets_Z   = targets_Z[condition]
	
	return targets_R, targets_Z, active

#==============================================================================
# FUNCTION - save_intersection
#==============================================================================
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
#==============================================================================
def save_intersection(targets_R, targets_Z, intersect_R, intersect_Z):

    intersect_R = np.append(intersect_R,targets_R)
    intersect_Z = np.append(intersect_Z,targets_Z)

    return intersect_R, intersect_Z

#==============================================================================
