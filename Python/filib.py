#==============================================================================
# MODULE - filib.py
#==============================================================================
#> 
#> DESCRIPTION:
#> Contains auxiliary functions to trace magnetic field lines
#>
#> REVISION HISTORY:
#> 16 03 2021 - Carlos Romero
#>
#> TO DO:
#> -
#> 
#==============================================================================
import numpy as np

from scipy.interpolate import interpn

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
# FUNCTION - step_phi_cylindrical_euler
#==============================================================================
#> 
#> DESCRIPTION:
#> Computes the next step of a finite differences equation by the Euler method.
#> 
#> INPUT/OUTPUT:
#> - [IN]    B_R_inter   - interpolated radial magnetic field for targets
#> - [IN]    B_phi_inter - interpolated toroidal magnetic field for targets
#> - [IN]    B_Z_inter   - interpolated vertical magnetic field for targets
#> - [IN]    delta_phi   - increment in toroidal direction
#> - [INOUT] targets_R   - R-coordinates of target points
#> - [INOUT] targets_Z   - Z-coordinates of target points
#> 
#> EXAMPLE:
#> - Example: unnecessary
#>
#> REVISION HISTORY:
#> 16 03 2021 - Carlos Romero
#> 
#> TO DO:
#> -
#> 
#==============================================================================
def step_phi_cylindrical_euler(B_R_inter, B_phi_inter, B_Z_inter, delta_phi, targets_R, targets_Z):
		
	# - the order is crucial: first update the Z coordinate!
	targets_Z = targets_Z + delta_phi * ( B_Z_inter / B_phi_inter ) * targets_R
	targets_R = targets_R + delta_phi * ( B_R_inter / B_phi_inter ) * targets_R
		
	return targets_R, targets_Z
	
#==============================================================================
# FUNCTION - step_phi_cylindrical_rk4
#==============================================================================
#> 
#> DESCRIPTION:
#> Computes the next step of a finite differences equation by the Runge-Kutta
#> 4th order method.
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
#> 16 03 2021 - Carlos Romero
#> 
#> TO DO:
#> -
#> 
#==============================================================================
def step_phi_cylindrical_rk4(grid_points, B_R, B_phi, B_Z, targets_R, targets_phi, targets_Z, delta_phi):
		
	coeff_rk4 = np.array([0, 0., 0.5, 0.5, 1.]) * delta_phi 
	
	kZ = np.zeros((5,len(targets_R)))
	kR = np.zeros((5,len(targets_R)))
	
	kZ[0,:] = np.copy(targets_Z)
	kR[0,:] = np.copy(targets_R)
	
	for k in range(1,5):
	
		targets     = np.vstack([targets_phi, kZ[0,:] + coeff_rk4[k] * kZ[k-1,:], \
											  kR[0,:] + coeff_rk4[k] * kR[k-1,:]])
		targets     = np.transpose(targets)	
	
		# - np.nan as fill_value may be problematic
		B_phi_inter = interpn(grid_points, B_phi, targets, method='linear', bounds_error=False, fill_value=np.nan)
		B_Z_inter   = interpn(grid_points, B_Z,   targets, method='linear', bounds_error=False, fill_value=np.nan)
		B_R_inter   = interpn(grid_points, B_R,   targets, method='linear', bounds_error=False, fill_value=np.nan)
		
		kZ[k,:] = ( B_Z_inter / B_phi_inter ) * ( kR[0,:] + coeff_rk4[k] * kZ[k-1,:] )
		kR[k,:] = ( B_R_inter / B_phi_inter ) * ( kR[0,:] + coeff_rk4[k] * kR[k-1,:] )
		
	
	targets_Z = targets_Z + 1 / 6. * delta_phi * ( kZ[1,:] + 2 * kZ[2,:] + 2 * kZ[3,:] + kZ[4,:] )
	targets_R = targets_R + 1 / 6. * delta_phi * ( kR[1,:] + 2 * kR[2,:] + 2 * kR[3,:] + kR[4,:] ) 
		
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
def delete_outer_particles(rmin, rmax, zmin, zmax, targets_R, targets_Z, active):

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
#> - Example:
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
#==============================================================================
