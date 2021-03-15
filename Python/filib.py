#==============================================================================
# MODULE - filib.py
#==============================================================================
#> 
#> DESCRIPTION:
#> Contains auxiliary functions to trace magnetic field lines
#>
#> REVISION HISTORY:
#> 15 03 2021 - Carlos Romero
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
#> 10 03 2021 - Carlos Romero
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

    if mode == 'line-vertical':
    
    	targets_R = (rmax + rmin) / 2
    	targets_Z = np.linspace(zmin,zmax,ntargets)
    	
    	print('\n')
    	print('=====================================================')
    	print('        DEFAULT LINE WITHIN PROVIDED MARGINS:        ')
    	print('            R = (rmax + rmin) / 2 = ' + str(targets_R))
    	print('=====================================================')
    	print('\n')
    	
    elif mode == 'line-horizontal':
    	targets_R = np.linspace(rmin,rmax,ntargets)
    	targets_Z = (zmax + zmin) / 2
    	
    	print('\n')
    	print('=====================================================')
    	print('        DEFAULT LINE WITHIN PROVIDED MARGINS:        ')
    	print('            Z = (zmax + zmin) / 2 = ' + str(targets_Z))
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
# FUNCTION - evolve_cylindrical
#==============================================================================
#> 
#> DESCRIPTION:
#> Traces discretized field lines in cylindrical geometry: field line equation.
#> CAREFUL with the conversion from toroidal to cylindrical geometry!
#> 
#> INPUT/OUTPUT:
#> - [IN]  array0      - particle's initital positions ([phi, Z, R])
#> - [IN]  step_phi_mm - trajectory's step size in mm
#> - [IN]  B_inter     - interpolated magnetic field for targets ([phi, Z, R])
#> - [OUT] array1      - particle's update position ([phi, Z, R])
#> 
#> EXAMPLE:
#> - Example: unnecessary
#>
#> REVISION HISTORY:
#> 15 03 2021 - Carlos Romero
#> 
#> TO DO:
#> -
#> 
#==============================================================================
def evolve_cylindrical(array0, step_phi_mm, B_inter):

	# - the array components are:
	# phi0 = array0[0]
	# Z0   = array0[1]
	# R0   = array0[2]
	
	# - analogously for the magnetic field components:
	# B_phi_inter = B_inter[0]
	# B_Z_inter   = B_inter[1]
	# B_R_inter   = B_inter[2]
	
	# - we compute the angle that a given particle (R0) must run to
	# have a given step size in mm
	delta_phi = step_phi_mm / ( array0[2] * 1000 )
	
	# - declare array to be updated
	array1 = np.empty(np.shape(array0))

	# - update the particle's position: phi coordinate is symmetric
	array1[0] = (array0[0] + delta_phi ) % ( 2 * np.pi )
	array1[1] = array0[1] + delta_phi * ( B_inter[1] / B_inter[0] ) * array0[2]
	array1[2] = array0[2] + delta_phi * ( B_inter[2] / B_inter[0] ) * array0[2]

	return array1

#==============================================================================
# FUNCTION - inside_boundaries
#==============================================================================
#> 
#> DESCRIPTION:
#> Computes boolean of particles position with respect to boundaries: 
#> inside (1) or outside (0).
#> 
#> INPUT/OUTPUT:
#> - [IN]  rmin        - minimum radial coordinate
#> - [IN]  rmax        - maximum radial coordinate
#> - [IN]  zmin        - minimum vertical coordinate
#> - [IN]  zmax        - maximum vertical coordinate
#> - [IN]  array1      - particle's update position ([phi, Z, R])
#> - [OUT] inside_bool - boolean with information regarding if particle is 
#>                     still inside (1) boundaries or outside (0)
#> 
#> EXAMPLE:
#> - Example: omitting unnecessary variables for illustration
#>      rmin = 2.8
#>      array1 = [2. -1.3 2.7]
#>      inside_bool = inside_boundaries(array1, ~, rmin, ~, ~, ~)
#>      $ inside_bool = 0
#>
#> REVISION HISTORY:
#> 15 03 2021 - Carlos Romero
#> 
#> TO DO:
#> - 
#> 
#==============================================================================
def inside_boundaries(rmin, rmax, zmin, zmax, array1):

	# - the array components are:
	# phi0 = array0[0]
	# Z0   = array0[1]
	# R0   = array0[2]

	condition_R = (rmin <= array1[2]) & (array1[2] <= rmax)	
	condition_Z = (zmin <= array1[1]) & (array1[1] <= zmax)
	
	condition   = condition_R & condition_Z
	
	inside_bool = condition * 1
	
	return inside_bool

#==============================================================================
# FUNCTION - inside_limiter
#==============================================================================
#> 
#> DESCRIPTION:
#> Deletes particles that are outside simulation boundaries.
#> 
#> INPUT/OUTPUT:
#> - [IN]  array1        - particle's update position ([phi, Z, R])
#> - [IN]  limiter_inter - interpolated limiter information of targets
#> - [OUT] inside_bool   - boolean with information regarding if particle is 
#>                       still inside (1) boundaries or outside (0)
#> 
#> EXAMPLE:
#> - Example:
#>
#> REVISION HISTORY:
#> 15 03 2021 - Carlos Romero
#> 
#> TO DO:
#> - 
#> 
#==============================================================================
def inside_limiter(limiter_inter, array1):

	# - if the particle is close to the wall, its interpolation value is
	# between 0 and 1
	condition = limiter_inter != 0.
	
	inside_bool = condition * 1
	
	return inside_bool

#==============================================================================
# FUNCTION - metric_cylindrical
#==============================================================================
#> 
#> DESCRIPTION:
#> Computes length element of trajectory in cylindrical coordinates.
#> 
#> INPUT/OUTPUT:
#> - [IN]  array0 - particle's original position ([phi, Z, R])
#> - [IN]  array1 - particle's update position ([phi, Z, R])
#> - [OUT] metric - differential length element 
#> 
#> EXAMPLE:
#> - Example:
#>
#> REVISION HISTORY:
#> 15 03 2021 - Carlos Romero
#> 
#> TO DO:
#> - 
#> 
#==============================================================================
def metric_cylindrical(array0, array1):

	# - the array components are:
	# phi0 = array0[0]
	# Z0   = array0[1]
	# R0   = array0[2]

	# - differentials and mean of radius between phi = k and phi = k-1
	delta_phi = array1[0] - array0[0]
	delta_Z   = array1[1] - array0[1]
	delta_R   = array1[2] - array0[2]
	mean_R  = (array1[2] + array0[2]) / 2 
	
	# - metric in cyllindrical coordinates between k and k-1 with element-wise
	# operations
	metric = np.sqrt( delta_R ** 2 + delta_Z ** 2 + mean_R ** 2 * delta_phi ** 2 )			
	
	return metric
	
#==============================================================================
#==============================================================================
