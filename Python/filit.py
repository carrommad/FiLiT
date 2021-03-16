#==============================================================================
# PROGRAM - filit.py
#==============================================================================
#> 
#> DESCRIPTION:
#> Contains necessary constants and subroutines to trace field lines in 2- or 
#> 3-dimensional systems to visualize flow structures.
#>
#> Tracing field lines requires:
#> - starting point of the curve: consider a line or a matrix of points
#> - vector field considered pointwise tangent to the curve
#> 
#> As we are dealing with a discretized space, both the space and the vector 
#> field are given in a multidimensional grid. To compute a realistic field
#> line, we will need to do:
#> - multidimensional interpolation to obtain precise field values
#> - coordinate transformations to account for specific geometries
#> 
#> REVISION HISTORY:
#> 16 03 2021 - Carlos Romero
#> 
#> TO DO:
#> - Comment 'plt.show()' for programmed shutdown
#> - In RK4, when doing interpolation: careful with 'fill_value=np.nan'!
#> - 
#> 
#==============================================================================


#==============================================================================
# 0.                            M O D U L E S
#==============================================================================
import netCDF4           as nc
import numpy             as np

import time

from matplotlib import pyplot as plt
from filib      import read_data_hint,                 \
				       read_limiter,                   \
				       limiter_boundary,               \
				       create_targets,                 \
				       step_phi_cylindrical_euler,     \
				       step_phi_cylindrical_rk4,       \
				       delete_outer_particles,         \
				       delete_outer_particles_limiter

from progress.bar import Bar

from scipy.interpolate import interpn


#==============================================================================
# 1.                             I N P U T S
#==============================================================================	
tar_mode     = 'matrix'#'line-horizontal' # Target's distribution: 'line-vertical/horizontal' - 'matrix'

ntargets1D   =  9                # Number of targets per line

nturns       =  50              # Number of turns to the reactor

step_phi_mm  =  80                # Toroidal step in mm, for major radius R0

filepath     = 'snapfile.80.nc'  # Full name of the input file (HINT output)

limpath      = 'limiter.nc'      # Full name of the limiter/divertor geometry file

lim_mode     =  1                # Simulation boundaries (0), real limiter (1) 


#==============================================================================
# 2.                           W E L C O M E
#==============================================================================	
print('\n')
print('=====================================================')
print('                                                     ')
print(' WELCOME TO                                          ')
print('                                                     ')
print('     .......... ...  ....       ... ............     ')
print('     :::::::::: :::  ::::       ::: ::::::::::::     ')
print('     ::::            ::::               ::::         ') 
print('     ::::...    ...  ::::       ...     ::::         ')
print('     :::::::    :::  ::::       :::     ::::         ')
print('     ::::       :::  ::::       :::     ::::         ')
print('     ::::       :::  :::::::::: :::     ::::         ')
print('     ····       ···  ·········· ···     ····   v1.1  ')
print('                                                     ')
print('=====================================================')

print('\n')
print('=====================================================')
print('                     I N P U T S                     ') 
print('=====================================================')
print('   - target mode                ->  ' + tar_mode      )
print('   - number of targets per line ->  ' + str(ntargets1D))
print('   - number of turns            ->  ' + str(nturns)    )
print('   - toroidal step size in mm   ->  ' + str(step_phi_mm))
print('   - full input file name       ->  ' + filepath       )
print('   - full limiter file name     ->  ' + limpath       )
print('   - limiter mode               ->  ' + str(lim_mode) )
print('=====================================================')


#==============================================================================
# 3.                      R E A D   D A T A   H I N T
#==============================================================================	
# - magnetic equilibrium dataset
dataset = nc.Dataset(filepath)

# - read and extract data
mtor, rminb, rmaxb, zminb, zmaxb, B_R_original, B_phi_original, B_Z_original = read_data_hint(dataset)

# - make copy of original data for possible replication
B_R   = np.copy(B_R_original)
B_phi = np.copy(B_phi_original)
B_Z   = np.copy(B_Z_original)

# - replicate symmetric toroidal data 'mtor - 1' number of times, to have
# a copy of mtor-times the original, completing the stellarator toroidally
for n in range(mtor-1):
	B_R   = np.append(B_R, B_R_original, axis = 0)
	B_phi = np.append(B_phi, B_phi_original, axis = 0)
	B_Z   = np.append(B_Z, B_Z_original, axis = 0)

# - close dataset
dataset.close()
	

#==============================================================================
# 4.                   R E A D   D A T A   L I M I T E R
#==============================================================================	
# - limiter conditional
if lim_mode == 1:

	# - limiter/divertor geometry data
	limset = nc.Dataset(limpath)
	
	# - read and extract data
	limiter_original = read_limiter(limset)
	
	# - make copy of original data for possible replication
	limiter = np.copy(limiter_original)
	
	# - replicate symmetric toroidal data 'mtor - 1' number of times, to have
	# a copy of mtor-times the original, completing the stellarator toroidally
	for n in range(mtor-1):
		limiter = np.append(limiter, limiter_original, axis = 0)
	
	# - close dataset
	limset.close()


#==============================================================================
# 5.         I N P U T   D I M E N S I O N S   A N D   G R I D
#==============================================================================	
# - extracting dimensions from data
grid_size  = B_R.shape

# - recreating original grid: careful with the orderof original grids!
# - original order [r,z,phi] but it appears transposed: [phi,z,r]
n_grid_phi = grid_size[0] 
n_grid_Z   = grid_size[1]
n_grid_R   = grid_size[2]


#==============================================================================
# 6.            G R I D   D E N S I T Y   F R O M   S T E P  (phi)
#==============================================================================	
# - computing major radius 
R0 = (rminb + rmaxb) / 2

# - for a given toroidal step, compute grid density only if bigger than
# original
n_grid_phi_enhanced = max(n_grid_phi, int(2 * np.pi * R0 * 1000 / step_phi_mm ))

# - conditional for a too big toroidal step
if n_grid_phi_enhanced == n_grid_phi:
	print('\n')
	print('=====================================================')
	print('      ATTENTION: USING ORIGINAL GRID SIZE: ' + str(n_grid_phi))
	print('=====================================================')
	print('\n')

# - printing grid density info if step size actually enhances original grid
else:
	print('\n')
	print('=====================================================')
	print('    INFO: ENHANCED GRID vs. ORIGINAL: ' + str(n_grid_phi_enhanced) + ' > ' + str(n_grid_phi))
	print('=====================================================')
	print('\n')


#==============================================================================
# 7.                  G R I D S   A N D   S T E P S  (R, Z)
#==============================================================================	
# - recreating original grid of periodic dimension
grid_phi = np.linspace(0,     2*np.pi, n_grid_phi+1)
grid_phi = grid_phi[0:n_grid_phi]

# - recreating original grids of non-periodic dimensions
grid_Z   = np.linspace(zminb, zmaxb,   n_grid_Z)
grid_R   = np.linspace(rminb, rmaxb,   n_grid_R)

# - computing steps of regular grids
step_Z = (zmaxb - zminb ) / (n_grid_Z - 1)
step_R = (rmaxb - rminb ) / (n_grid_R - 1)

# - we provide a different name to the toroidal step as it is also a
# differential for the connection length computation and to avoid confusion 
# with the imported function 'step_phi_cylindrical_xxxxx'
# - due to periodicity, the grid size is differently defined
delta_phi = 2 * np.pi / n_grid_phi_enhanced


#==============================================================================
# 8.                      I N I T I A L   T A R G E T S
#==============================================================================
# - define a radial left and right margins to create targets far from
# boundaries
#left  = rminb + 140 * step_R
#right = rmaxb - 200 * step_R
print('\n')
print('NOTE: ad hoc target horizontal positions! ')
left  = 3.23
right = 3.74

# - analagous for vertical margins
bottom = zminb + 70 * step_Z
top    = zmaxb - 70 * step_Z

# - creating targets for given input values
targets_R, targets_Z = create_targets(tar_mode, ntargets1D, left, right, bottom, top)

# - depending on the selected target mode we have ntargets1D or ntargets1D^2
# total targets
ntargets = len(targets_R)


#==============================================================================
# 9.   N U M B E R   O F   S T E P S   A N D   I N T E R S E C T I O N
#==============================================================================
# - iterate phi a certain number of steps given by grid and toroidal turn
nsteps = nturns * n_grid_phi_enhanced

# - declare array to save intersections
print('\n')
print('NOTE: saving intersections on symmetric toroidal planes! ')
# NOTE: multiplying by mtor allows us to save intersection points at every symmetry toroidal angle
intersect_R = np.nan * np.empty((ntargets, nturns * mtor))
intersect_Z = np.nan * np.empty((ntargets, nturns * mtor))

# - assing initial targets to 0-th column of intersections
intersect_R[:, 0] = np.copy(targets_R)
intersect_Z[:, 0] = np.copy(targets_Z)


#==============================================================================
# 10.           F O R   3-D   I N T E R P O L A T I O N   G R I D
#==============================================================================
# - create array of grid points in space
grid_points = (grid_phi, grid_Z, grid_R)


#==============================================================================
# 11.              F O R   C O N N E C T I O N   L E N G T H
#==============================================================================
# - initial definition of connection length 
connection_length = np.zeros(np.shape(targets_R))

# - full copy of targets in the present iteration step in the 'for' loop
targets_R_lc1 = np.copy(targets_R)
targets_Z_lc1 = np.copy(targets_Z)


#==============================================================================
# 12.              F O R   P A R T I C L E   T R A C K I N G
#==============================================================================
# - full copy of targets at previous iteration step in the 'for' loop
targets_R_lc0 = np.copy(targets_R)
targets_Z_lc0 = np.copy(targets_Z)

# - define all particles initially as 'active' = 1
active = np.ones(np.shape(targets_R))

# - plane counter for intersection saving
plane = 1


#==============================================================================
# 13.                   M A I N   L O O P :   T U R N S
#==============================================================================
# - process bar for loop
print('\n')
bar = Bar('Processing turns', max=nturns)

# - toroidal loop
for step in range(nsteps):
		
				
	#==============================================================================
	# a.            F U L L   T A R G E T S   P R E V I O U S   S T E P
	#==============================================================================
	# - full copy of particles of previous iteration, modifying only active
	# particles
	targets_R_lc0[active == 1.] = np.copy(targets_R)
	targets_Z_lc0[active == 1.] = np.copy(targets_Z)
		
		
	#==============================================================================
	# b.             T O R O I D A L   S T E P   A N D   A N G L E
	#==============================================================================
	# - fraction of the angle 
	fraction = step % n_grid_phi_enhanced
	
	# - angle in radians modulo 2*pi
	angle = fraction * delta_phi
	
	
	#==============================================================================
	# c.    D E F I N E   T A R G E T S   F O R   I N T E R P O L A T I O N
	#==============================================================================
	# - create array of phi position of targets at this step
	targets_phi = angle * np.ones(np.shape(targets_Z))
	
	"""
	# - stack together [phi,Z,R] to obtain 3D positions
	targets = np.vstack([targets_phi, targets_Z, targets_R])
	
	# - transpose to coincide with dimensions of fields
	targets = np.transpose(targets)		
	

	#==============================================================================
	# d.         3-D   I N T E R P O L A T I O N   G R I D :   F I E L D S
	#==============================================================================
	# - interpolate magnetic fields in 3D (it could also be done with the manually
	# defined interpolation functions (2D and 3D) in 'field_line_mod'
	B_phi_inter = interpn(grid_points, B_phi, targets, method='linear')
	B_Z_inter   = interpn(grid_points, B_Z,   targets, method='linear')
	B_R_inter   = interpn(grid_points, B_R,   targets, method='linear')
	
	
	#==============================================================================
	# e.          T O R O I D A L   F I N I T E   D I F F E R E N C E S
	#==============================================================================
	# - evolve finite differences with Euler of 1st order
	targets_R, targets_Z = step_phi_cylindrical_euler(B_R_inter, B_phi_inter, B_Z_inter, delta_phi, targets_R, targets_Z)
	"""
	
	#==============================================================================
	# ?.             T E S T I N G :   R U N G E - K U T T A    4  
	#==============================================================================
	# - evolve finite differences with Runge-Kutta of 4th order
	targets_R, targets_Z = step_phi_cylindrical_rk4(grid_points, B_R, B_phi, B_Z, targets_R, targets_phi, targets_Z, delta_phi)

	
	#==============================================================================
	# f.   D E L E T E   P A R T I C L E S   O U T S I D E   B O U N D A R I E S
	#==============================================================================
	# - deleting particles outside numerical boundaries 
	targets_R, targets_Z, active = delete_outer_particles(rminb, rmaxb, zminb, zmaxb, targets_R, targets_Z, active)
	
	
	#==============================================================================
	# g.     D E L E T E   P A R T I C L E S   O U T S I D E   L I M I T E R
	#==============================================================================
	# - conditional limiter mode
	if lim_mode == 1:
	
	
		#==============================================================================
		# *     R E D E F I N E   T A R G E T S   F O R   I N T E R P O L A T I O N
		#==============================================================================
		# - create array of phi position of targets at this step
		targets_phi = angle * np.ones(np.shape(targets_Z))
	
		# - stack together [phi,Z,R] to obtain 3D positions and transpose to coincide
		# with dimensions of fields
		targets = np.vstack([targets_phi, targets_Z, targets_R])
		targets = np.transpose(targets)
		
		
		#==============================================================================
		# *        3-D   I N T E R P O L A T I O N   G R I D :   L I M I T E R
		#==============================================================================
		# - interpolate limiter
		limiter_inter = interpn(grid_points, limiter, targets)#, method='linear', bounds_error=True, fill_value=0.)
	
		# - delete particles outside limiter/divertor geometry
		targets_R, targets_Z, active = delete_outer_particles_limiter(limiter_inter, targets_R, targets_Z, active)
	
	
	#==============================================================================
	# h.   C O M P U T E   D I S T A N C E S   B E T W E E N   P H I   S T E P S
	#==============================================================================
	# - full copy of only active targets prior to deletion for connetion length
	# computation
	targets_R_lc1[active == 1.] = np.copy(targets_R)
	targets_Z_lc1[active == 1.] = np.copy(targets_Z)

	# - differentials and mean of radius between phi = k and phi = k-1
	delta_R =  targets_R_lc1 - targets_R_lc0
	delta_Z =  targets_Z_lc1 - targets_Z_lc0
	mean_R  = (targets_R_lc1 + targets_R_lc0) / 2 
	
	# - metric in cyllindrical coordinates between k and k-1 with element-wise
	# operations
	metric = np.sqrt( delta_R ** 2 + delta_Z ** 2 + mean_R ** 2 * delta_phi ** 2 )			

	# - summing to prior connection length for active particles only
	connection_length = connection_length + metric * active
	
	
	#==============================================================================
	# i.                 E X I T   I F   N O   T A R G E T S
	#==============================================================================
	# - exit loop if there are no more targets
	if len(targets_R) == 0:
		print('\n')
		print('=====================================================')	
		print('        ATTENTION: YOU RUN OUT OF PARTICLES!         ')
		print('=====================================================')
		break
		
		
	#==============================================================================
	# j.            S A V I N G   I N T E R S E C T I O N S
	#==============================================================================
	# - commented option to use mtor-symmetry of trajectories to save more 
	# intersection points: more efficient (you would need to indent 'if' and 
	# 'bar.next')
	if (fraction % (n_grid_phi_enhanced / mtor) == 0 and step > 0):
	
		# - save positions
		intersect_Z[:, plane] = np.copy(targets_Z_lc1)
		intersect_R[:, plane] = np.copy(targets_R_lc1)
		
		# - update intersection plane counter
		plane = plane + 1
		
		# - if full turn completed (indent/unindent)
		if (fraction == 0 and step > 0):

			# - step in process bar (indent/unindent)
			bar.next()
		

#==============================================================================
# 14.                E N D   M A I N   L O O P :   T U R N S
#==============================================================================	
# - step in process bar
bar.next()

# - ending bar counter and printing
bar.finish()


#==============================================================================
# 15.                 D A T A   S A V I N G   P A T H 
#==============================================================================
# - folder name
folder = 'results'

# - full string for saving data
str_save = folder + '/' + tar_mode + '_ntargets-' + str(ntargets) + '_nturns-' + str(nturns)


#==============================================================================
# 16.                S A V I N G   I N T E R S E C T I O N S   
#==============================================================================
# - create a netCDF dataset to save data
saved_data = nc.Dataset(str_save + '.nc', 'w')

# - declare number of particles and turns
saved_data.createDimension('Np', ntargets)
saved_data.createDimension('Nt', nturns*mtor)


# - declare variables
plot_intersect_R = saved_data.createVariable('plot_intersect_R',  'f4', ('Np','Nt',))#,  zlib=True)
plot_intersect_Z = saved_data.createVariable('plot_intersect_Z',  'f4', ('Np','Nt',))#,  zlib=True)

# - assign values to variables
plot_intersect_R[:] = intersect_R
plot_intersect_Z[:] = intersect_Z


#==============================================================================
# 17.        S A V I N G   C O N N E C T I O N   L E N G T H   D A T A   
#==============================================================================
# - declare dimensions
saved_data.createDimension('R', ntargets)

# - declare variables
plot_connection_length = saved_data.createVariable('plot_connection_length',  'f4', ('Np',))#,  zlib=True)

# - assign values to netcdf variables
plot_connection_length[:] = connection_length


#==============================================================================
# 18.         S A V I N G   L I M I T E R   B O U N D A R I E S
#==============================================================================
# - conditional limiter mode
if lim_mode == 1:

	# - consider outer points as 1
	limiter_R, limiter_Z = limiter_boundary(grid_R, grid_Z, limiter[0,:,:])

	# - declare number of intersections
	saved_data.createDimension('L', len(limiter_R))
	
	# - declare variables
	plot_limiter_R = saved_data.createVariable('plot_limiter_R',  'f4', ('L',))#,  zlib=True)
	plot_limiter_Z = saved_data.createVariable('plot_limiter_Z',  'f4', ('L',))#,  zlib=True)

	# - assign values to variables
	plot_limiter_R[:] = limiter_R
	plot_limiter_Z[:] = limiter_Z

# - close netCDF data to be saved
saved_data.close()


#==============================================================================
# 19.               P L O T   1 :   I N T E R S E C T I O N S
#==============================================================================
# - plotting intersections
plt.figure(1)

# - plotting loop for intersections
for target in range(ntargets):

	plt.plot(intersect_R[target,:], intersect_Z[target,:], 'o', markersize=0.5)
	
# - plotting initial targets
plt.plot(intersect_R[:,0], intersect_Z[:,0], 'g*', markersize=1.5)
plt.title('Field line intersections at $\phi$ = 0')
plt.xlabel('R (m)')
plt.ylabel('Z (m)')
plt.grid(True)


#==============================================================================
# 20.                  P L O T   1 :   L I M I T E R
#==============================================================================
# - conditional limiter mode
if lim_mode == 1:

	# - plot limiter line
	plt.plot(limiter_R, limiter_Z, 'ko', markersize=0.5)

# - save figure
plt.savefig(str_save + '.png', dpi=1000)

# - show figure; if the simulation is programmed: comment 'plt.show()' and 
# uncomment 'plt.close()'
plt.show() 
#plt.close()


#==============================================================================
# 21.           P L O T   2 :   C O N N E C T I O N   L E N G T H
#==============================================================================
# - plotting connection length in logarithmic scale
plt.figure(2)

# - print lines with different colors
if tar_mode == 'matrix':

	for line in range(int(ntargets / ntargets1D)):

		# - indices for plotting data
		start_index = line * ntargets1D
		end_index   = line * ntargets1D + ntargets1D
		
		plt.plot(intersect_R[start_index:end_index,0], connection_length[start_index:end_index], '-')
		plt.xlabel('R (m)')
		plt.title('Radial connection length')

# - for vertical line
elif 'horizontal' in tar_mode:

	plt.plot(intersect_R[:,0], connection_length, '-')
	plt.xlabel('R (m)')
	plt.title('Radial connection length')
		
# - for horizontal line
elif 'vertical' in tar_mode:

	plt.plot(intersect_Z[:,0], connection_length, '-')
	plt.xlabel('Z (m)')
	plt.title('Vertical connection length')

plt.yscale('log')
plt.ylabel('L_c (m)')

# - save figure
plt.savefig(str_save + '_Lc.png', dpi=1000)

# - show figure; if the simulation is programmed: comment 'plt.show()' and 
# uncomment 'plt.close()'
plt.show() 
#plt.close()


#==============================================================================
# 22.                          G O O D B Y E
#==============================================================================
print('\n')
print('=====================================================')
print('            FIELD LINE TRACING FINISHED!             ')
print('=====================================================')
print('\n')


#==============================================================================
#==============================================================================

