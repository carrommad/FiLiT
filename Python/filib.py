#==============================================================================
# MODULE - filib.py
#==============================================================================
#> 
#> DESCRIPTION:
#> Contains auxiliary functions to trace magnetic field lines
#>
#> REVISION HISTORY:
#> 17 03 2021 - Carlos Romero
#>
#> TO DO:
#> -
#> 
#==============================================================================
import numpy as np

from scipy.interpolate import interpn

from matplotlib        import pyplot as plt


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
#> - [OUT] B_phi   - 3-tensor of toroidal magnetic field in the reactor
#> - [OUT] B_Z     - 3-tensor of vertical magnetic field in the reactor
#> - [OUT] B_R     - 3-tensor of radial magnetic field in the reactor
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

    return mtor, rminb, rmaxb, zminb, zmaxb, B_phi, B_Z, B_R
    
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
#> - [IN]  grid_Z        - Z-coordinates of grid points
#> - [IN]  grid_R        - R-coordinates of grid points
#> - [IN]  limiter_slice - limiter toroidal slice 
#> - [OUT] limiter_Z     - Z-coordinates of limiter line points
#> - [OUT] limiter_R     - R-coordinates of limiter line points
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
def limiter_boundary(grid_Z, grid_R, limiter_slice):
	
	# - compute slice dimensions
	Nz = limiter_slice.shape[0]
	Nr = limiter_slice.shape[1]
	
	# - declare lists to be appended
	limiter_Z = []
	limiter_R = []
	
	# - looping limiter values 
	for i in range(Nz):
		for j in range(Nr-1):
					
			# - stepping inside the limiter
			if limiter_slice[i,j+1] - limiter_slice[i,j] == 1:
			
				limiter_Z.append(grid_Z[i])
				limiter_R.append(grid_R[j])
				
			# - stepping outside the limiter
			elif limiter_slice[i,j+1] - limiter_slice[i,j] == -1:
				
				limiter_Z.append(grid_Z[i])
				limiter_R.append(grid_R[j+1])
			
	# - creating numpy arrays	
	limiter_Z = np.array(limiter_Z)
	limiter_R = np.array(limiter_R)
	
	return limiter_Z, limiter_R

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
#> - [OUT] targets_Z - Z-coordinates of target particles
#> - [OUT] targets_R - R-coordinates of target particles
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
    
    	targets_Z = np.linspace(zmin,zmax,ntargets)
    	targets_R = (rmax + rmin) / 2
    	
    	print('\n')
    	print('=====================================================')
    	print('        DEFAULT LINE WITHIN PROVIDED MARGINS:        ')
    	print('            R = (rmax + rmin) / 2 = ' + str(targets_R))
    	print('=====================================================')
    	
    elif mode == 'line-horizontal':
    
    	targets_Z = (zmax + zmin) / 2
    	targets_R = np.linspace(rmin,rmax,ntargets)
    	
    	print('\n')
    	print('=====================================================')
    	print('        DEFAULT LINE WITHIN PROVIDED MARGINS:        ')
    	print('            Z = (zmax + zmin) / 2 = ' + str(targets_Z))
    	print('=====================================================')
    	
    elif mode == 'matrix':
    	targets_Z = np.linspace(zmin,zmax,ntargets)
    	targets_R = np.linspace(rmin,rmax,ntargets)
    	
    targets_R, targets_Z = np.meshgrid(targets_R, targets_Z)
    
    targets_Z = targets_Z.flatten()
    targets_R = targets_R.flatten()
    
    return targets_Z, targets_R
	
#==============================================================================
# FUNCTION - step_phi_cylindrical_rk4
#==============================================================================
#> 
#> DESCRIPTION:
#> Computes the next step of a finite differences equation by the Runge-Kutta
#> 4th order method.
#> 
#> INPUT/OUTPUT:
#> - [IN]    grid_points - original 3D grid for interpolation
#> - [IN]    B_phi       - toroidal magnetic field for targets
#> - [IN]    B_Z         - vertical magnetic field for targets
#> - [IN]    B_R         - radial magnetic field for targets
#> - [IN]    targets_phi - phi-coordinates of target points
#> - [INOUT] targets_Z   - Z-coordinates of target points
#> - [INOUT] targets_R   - R-coordinates of target points
#> - [IN]    delta_phi   - toroidal increment
#> 
#> EXAMPLE:
#> - Example: unnecessary
#>
#> REVISION HISTORY:
#> 17 03 2021 - Carlos Romero
#> 
#> TO DO:
#> -
#> 
#==============================================================================
def step_phi_cylindrical_rk4(grid_points, B_phi, B_Z, B_R, targets_phi, targets_Z, targets_R, delta_phi):
		
	# - Runge-Kutta coefficients: the 0th element has been added for completion
	# of the loop
	coeff_rk4 = np.array([0, 0., 0.5, 0.5, 1.]) * delta_phi 
	
	# - Runge-Kutta iterations: the 0th element has been added for completion
	# of the loop
	kZ = np.zeros((5,len(targets_R)))
	kR = np.zeros((5,len(targets_R)))
		
	# - loop to obtain Runge-Kutta evaluations
	for k in range(1,5):
		
		# - define Runge-Kutta values for Z and R 
		value_Z = coeff_rk4[k] * kZ[k-1,:]
		value_R = coeff_rk4[k] * kR[k-1,:] 
	
		# - evaluate function at value given by Runge-Kutta
		kZ[k,:], kR[k,:] = field_line_function(grid_points, B_phi, B_Z, B_R, targets_phi, targets_Z, targets_R, value_Z, value_R)

		
		
	# - Runge-Kutta finite differences of 4th order 
	targets_Z = targets_Z + 1 / 6. * delta_phi * ( kZ[1,:] + 2 * kZ[2,:] + 2 * kZ[3,:] + kZ[4,:] )
	targets_R = targets_R + 1 / 6. * delta_phi * ( kR[1,:] + 2 * kR[2,:] + 2 * kR[3,:] + kR[4,:] ) 
		
	return targets_Z, targets_R
	
#==============================================================================
# FUNCTION - delete_outer_particles
#==============================================================================
#> 
#> DESCRIPTION:
#> Deletes particles that are outside simulation boundaries.
#> 
#> INPUT/OUTPUT:
#> - [IN]  grid_points - original 3D grid for interpolation
#> - [IN]  B_phi       - toroidal magnetic field for targets
#> - [IN]  B_Z         - vertical magnetic field for targets
#> - [IN]  B_R         - radial magnetic field for targets
#> - [IN]  targets_phi - phi-coordinates of target points
#> - [IN]  targets_Z   - Z-coordinates of target points
#> - [IN]  targets_R   - R-coordinates of target points
#> - [IN]  value_Z     - vertical position to be evaluated
#> - [IN]  value_R     - radial position to be evaluated
#> - [OUT] output_Z    - result of evaluated field line computation in Z
#> - [OUT] output_R    - result of evaluated field line computation in R
#> 
#> EXAMPLE:
#> - Example: 
#>
#> REVISION HISTORY:
#> 17 03 2021 - Carlos Romero
#> 
#> TO DO:
#> - 
#> 
#==============================================================================
def field_line_function(grid_points, B_phi, B_Z, B_R, targets_phi, targets_Z, targets_R, value_Z, value_R):

	# - see Runge-Kutta definition
	targets = np.vstack([targets_phi, targets_Z + value_Z, targets_R + value_R])
	targets = np.transpose(targets)	
	
	# - interpolate fields: np.nan as fill_value may be problematic
	B_phi_inter = interpn(grid_points, B_phi, targets, method='linear', bounds_error=False, fill_value=1e-10)#np.nan)
	B_Z_inter   = interpn(grid_points, B_Z,   targets, method='linear', bounds_error=False, fill_value=1.)#np.nan)
	B_R_inter   = interpn(grid_points, B_R,   targets, method='linear', bounds_error=False, fill_value=1.)#np.nan)
	 
	# - compute function:
	output_Z = ( B_Z_inter / B_phi_inter ) * (targets_R + value_Z)
	output_R = ( B_R_inter / B_phi_inter ) * (targets_R + value_R)

	return output_Z, output_R
	
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
#> - [INOUT] targets_Z - Z-coordinates of target points
#> - [INOUT] targets_R - R-coordinates of target points
#> - [INOUT] active    - array containing information (0/1) regarding if 
#>                     particle is still inside boundaries
#> 
#> EXAMPLE:
#> - Example: omitting unnecessary variables for illustration
#>      rmin = 0
#>      targets_R = [2. 8.3 -1.3 2.7]
#>      targets_R = delete_outer_particles(rmin, ~, ~, ~, ~, targets_R, ~)
#>      $ targets_R = [2. 8.3 2.7]
#>
#> REVISION HISTORY:
#> 24 02 2021 - Carlos Romero
#> 
#> TO DO:
#> - 
#> 
#==============================================================================
def delete_outer_particles(rmin, rmax, zmin, zmax, targets_Z, targets_R, active):

	condition_Z = (zmin <= targets_Z) & (targets_Z <= zmax)
	condition_R = (rmin <= targets_R) & (targets_R <= rmax)	
	
	condition   = condition_R & condition_Z
	
	# - changing values of active particles only
	active[active == 1.] = condition
	
	targets_Z   = targets_Z[condition]
	targets_R   = targets_R[condition]
	
	return targets_Z, targets_R, active

#==============================================================================
# FUNCTION - delete_outer_particles_limiter
#==============================================================================
#> 
#> DESCRIPTION:
#> Deletes particles that are outside simulation boundaries.
#> 
#> INPUT/OUTPUT:
#> - [IN]    limiter_inter - interpolated limiter information of targets
#> - [INOUT] targets_Z     - Z-coordinates of target points
#> - [INOUT] targets_R     - R-coordinates of target points
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
def delete_outer_particles_limiter(limiter_inter, targets_Z, targets_R, active):

	# - if the particle is close to the wall, its interpolation value is
	# between 0 and 1, so we take the particles that are clearly not outside (0)
	condition = limiter_inter != 0.	
	
	# - changing values of active particles only
	active[active == 1.] = condition
	
	targets_Z = targets_Z[condition]
	targets_R = targets_R[condition]
	
	return targets_Z, targets_R, active
	
#==============================================================================
# FUNCTION - plot_Poincare_Lc
#==============================================================================
#> 
#> DESCRIPTION:
#> 
#> 
#> INPUT/OUTPUT:
#> - [IN] filename          - name with information on target mode (line/matrix)
#> - [IN] intersect_Z       - intersection with \phi=0 vertical positions, tensor
#>                          with rows = particles and columns = number of 
#>                          intersections
#> - [IN] intersect_R       - intersection with \phi=0 radial positions, tensor
#>                          with rows = particles and columns = number of 
#>                          intersections
#> - [IN] connection_length - connection length per particle
#> - [IN] limiter_Z         - verticla position of limiter transition (optional)
#> - [IN] limiter_R         - radial position of limiter transition (optional)
#> - [IN] plot_initial      - boolean for plotting initial target distribution
#> 
#> EXAMPLE:
#> - Example:
#>
#> REVISION HISTORY:
#> 17 03 2021 - Carlos Romero
#> 
#> TO DO:
#> - Programmed simulations: comment 'plt.show()', uncomment 'plt.close()'
#> - 
#> 
#==============================================================================
def plot_Poincare_Lc(filename, intersect_Z, intersect_R, connection_length, limiter_Z, limiter_R, plot_initial):

	#==============================================================================
	# a.                 C H E C K   T A R G E T   M O D E
	#==============================================================================
	# - number of target intersections
	ntargets = np.size(intersect_R[:,0])

	# - a way of seing if data is in a line or a matrix and computing number of 
	# targets per line
	if 'matrix' in filename:
		ntargets1D = int(np.sqrt(ntargets))
	elif 'line' in filename:
		ntargets1D = ntargets


	#==============================================================================
	# b.               P L O T   1 :   I N T E R S E C T I O N S
	#==============================================================================
	# - select figure
	plt.figure(1)

	# - plotting loop for targets intersecting with \phi = 0
	for target in range(1, ntargets):
		
		# - choosing divergin colormap
		cm = plt.cm.get_cmap('RdYlBu')
		
		# - color array depending on connection length value
		color = connection_length[target]*np.ones(np.shape(intersect_R[target,:]))
		
		# - marker size
		size = 0.2
		
		# - martker transparency
		transparency = 0.5
		
		# - colorbar limits
		low = min(connection_length)
		high = max(connection_length)
		
		# - plotting scattered points with color scale depending on connection length
		plt.scatter(intersect_R[target,:], intersect_Z[target,:], c=color, vmin=low, vmax=high, cmap=cm, s=size, alpha=transparency, marker='.')
		
		
	# - colorbar
	cbar = plt.colorbar()
	cbar.set_label("Lc (m)")

	# - conditional for plotting initial target distribution
	if plot_initial == 1:
		# - plotting initial targets
		plt.plot(intersect_R[:,0], intersect_Z[:,0], 'go', markersize=1.5)
	
	# - title, labels and grid
	plt.title('Field line intersections at $\phi$ = 0')
	plt.xlabel('R (m)')
	plt.ylabel('Z (m)')
	plt.grid(True)


	#==============================================================================
	# c.                    P L O T   1 :   L I M I T E R
	#==============================================================================
	# - plot limiter line
	plt.plot(limiter_R, limiter_Z, 'ko', markersize=1.5)


	#==============================================================================
	# d.                P L O T   1 :   S A V E   F I G U R E
	#==============================================================================
	# - save Poincaré plot
	plt.savefig(filename + '_Pp.png', dpi=1000)

	# - show/close figure
	plt.show() #plt.close()


	#==============================================================================
	# e.           P L O T   2 :   C O N N E C T I O N   L E N G T H
	#==============================================================================
	# - plotting connection length in logarithmic scale
	plt.figure(2)

	# - print lines with different colors for matrix target mode
	if 'matrix' in filename:

		# - loop over horizontal lines
		for line in range(int(ntargets / ntargets1D)):

			# - indices for plotting data
			start_index = line * ntargets1D
			end_index   = line * ntargets1D + ntargets1D
			
			# - radial intersection vs. connection length
			plt.plot(intersect_R[start_index:end_index,0], connection_length[start_index:end_index], '-')
			plt.xlabel('R (m)')
			plt.title('Radial connection length')

	# - for horizontal line
	elif 'horizontal' in filename:

		plt.plot(intersect_R[:,0], connection_length, '-')
		plt.xlabel('R (m)')
		plt.title('Radial connection length')
			
	# - for vertical line
	elif 'vertical' in tar_mode:

		plt.plot(intersect_Z[:,0], connection_length, '-')
		plt.xlabel('Z (m)')
		plt.title('Vertical connection length')

	# - common features
	plt.yscale('log')
	plt.ylabel('Lc (m)')

	# - save connection length plot
	plt.savefig(filename + '_Lc.png', dpi=1000)

	# - show/close figure
	plt.show() #plt.close()

#==============================================================================
#==============================================================================

