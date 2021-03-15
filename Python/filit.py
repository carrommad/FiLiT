#==============================================================================
# PROGRAM - field_line_tracer.py
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
#> 10 03 2021 - Carlos Romero
#> 
#> TO DO:
#> -
#> 
#==============================================================================


#==============================================================================
# 0.                            M O D U L E S
#==============================================================================
import netCDF4           as nc
import numpy             as np

import timeit
import time

from matplotlib     import pyplot as plt
from field_line_mod import read_data_hint,                 \
						   read_limiter,                   \
						   limiter_boundary,               \
						   create_targets,                 \
						   grid_neighbours_2D,             \
						   grid_neighbours_3D,             \
						   interpolate_field_2D,           \
						   interpolate_field_3D,           \
						   step_phi_cylindrical,           \
						   delete_outer_particles,         \
						   delete_outer_particles_limiter, \
						   save_intersection

from progress.bar import Bar

from scipy.interpolate import interpn, interp2d


#==============================================================================
# 1.                             I N P U T S
#==============================================================================	
tar_mode     = 'line'            # Target's distribution: 'line' - 'matrix'

ntargets1D   =  20               # Number of targets per line

nturns       =  300              # Number of turns to the reactor

step_phi_mm  =  1                # Toroidal step in mm, for major radius R0

filepath     = 'snapfile.80.nc'  # Full name of the input file (HINT output)

limpath      = 'limiter.nc'      # Full name of the limiter/divertor geometry file

lim_mode     =  1                # Simulation boundaries (0), real limiter (1) 

lc_mode      =  1                # Connection length computation for 'line': off (0), on (1)

debug        =  0                # Boolean for debug mode


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
print('     ····       ···  ·········· ···     ····   v1.0  ')
print('                                                     ')
print('=====================================================')

print('\n')
print('=====================================================')
print('                     I N P U T S                     ') 
print('=====================================================')
print('   - target mode                ->  ' + tar_mode           )
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
# a copy of mtor-times the original
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
	# a copy of mtor-times the original
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
# - defining original grids
grid_phi = np.linspace(0,     2*np.pi, n_grid_phi)
grid_Z   = np.linspace(zminb, zmaxb,   n_grid_Z)
grid_R   = np.linspace(rminb, rmaxb,   n_grid_R)

# - computing steps of regular grids
step_Z = (zmaxb - zminb ) / (n_grid_Z - 1)
step_R = (rmaxb - rminb ) / (n_grid_R - 1)

# - we provide a different name to the toroidal step as it is also a
# differential for the connection length computation and to avoid confusion 
# with the imported function 'step_phi'
delta_phi = 2 * np.pi / (n_grid_phi_enhanced - 1) 


#==============================================================================
# 8.                      I N I T I A L   T A R G E T S
#==============================================================================
# - define a margin to create targets far from boundaries
margin = 70 * step_R
	
# - creating targets for given input values
targets_R, targets_Z = create_targets(tar_mode, ntargets1D, rminb + margin, \
														rmaxb - margin, \
														zminb + margin, \
														zmaxb - margin)

# - depending on the selected target mode we have ntargets1D or ntargets1D^2
# total targets
ntargets = len(targets_R)


#==============================================================================
# 9.   N U M B E R   O F   S T E P S   A N D   I N T E R S E C T I O N
#==============================================================================
# - iterate phi a certain number of steps given by grid and toroidal turn
nsteps = nturns * n_grid_phi_enhanced

# - array to save intersections: initial values
intersect_R = np.copy(targets_R)
intersect_Z = np.copy(targets_Z)

# - for plotting
plot_targets_R = np.copy(targets_R)
plot_targets_Z = np.copy(targets_Z)


#==============================================================================
# 10.           F O R   3-D   I N T E R P O L A T I O N   G R I D
#==============================================================================
# - create array of grid points in space
grid_points = (grid_phi, grid_Z, grid_R)


#==============================================================================
# 11.              F O R   C O N N E C T I O N   L E N G T H
#==============================================================================
# - connection length conditional
if (lc_mode == 1 and tar_mode == 'line'):

	# - initial definition of connection length 
	connection_length = np.zeros(np.shape(targets_R))

	# - define all particles initially as 'active' = 1
	active = np.ones(np.shape(targets_R))

	# - full copy of targets at previous iteration step in the 'for' loop
	targets_R_lc0 = np.copy(targets_R)
	targets_Z_lc0 = np.copy(targets_Z)

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


#==============================================================================
# 13.                   M A I N   L O O P :   T U R N S
#==============================================================================
# - process bar for loop
bar = Bar('Processing turns', max=nturns)

# - toroidal loop
for step in range(nsteps):
		
				
	#==============================================================================
	# a.            F U L L   T A R G E T S   P R E V I O U S   S T E P
	#==============================================================================
	# - full copy of particles of previous iteration, modifying only active
	# particles
	targets_R_lc0[active == 1.] = targets_R
	targets_Z_lc0[active == 1.] = targets_Z
		
		
	#==============================================================================
	# b.             T O R O I D A L   S T E P   A N D   A N G L E
	#==============================================================================
	# - fraction of the angle 
	fraction = step % n_grid_phi_enhanced
	
	# - angle in radians modulo 2*pi
	angle = fraction * delta_phi
	
	
	#==============================================================================
	# ?.             D E B U G :   P L O T   T R A J E C T O R Y   
	#==============================================================================
	# - print every n_print iterations (roughly 10 points per poloidal turn)
	n_print = int(n_grid_phi_enhanced / mtor / 10)
	
	# - optional dynamic trajectory plot for debug
	if (debug == 1 and step % n_print == 0):
		plt.figure(3)
		if step == 0:
			plt.plot(targets_R_lc0, targets_Z_lc0,'bo')
		else:
			plt.plot(targets_R_lc0, targets_Z_lc0,'r.')
		
		plt.title('Iteration ' + str(step) + '/' + str(nsteps) + ' - Angle ' + str(int(100*angle)/100) + ' rad')
		plt.xlabel('R (m)')
		plt.ylabel('Z (m)')
		plt.axis([rminb, rmaxb, zminb, zmaxb])
		plt.pause(0.01)
			
	
	#==============================================================================
	# c.    D E F I N E   T A R G E T S   F O R   I N T E R P O L A T I O N
	#==============================================================================
	# - create array of phi position of targets at this step
	targets_phi = angle * np.ones(np.shape(targets_Z))
	
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
	# - evolve finite differences
	targets_R, targets_Z = step_phi_cylindrical(B_R_inter, B_phi_inter, B_Z_inter, delta_phi, rminb, rmaxb, targets_R, targets_Z)
	
	
	#==============================================================================
	# f.            F U L L   T A R G E T S   C U R R E N T   S T E P
	#==============================================================================
	# - connection length conditional
	if (lc_mode == 1 and tar_mode == 'line'):
	
		# - full copy of only active targets prior to deletion for connetion length
		# computation
		targets_R_lc1[active == 1.] = targets_R
		targets_Z_lc1[active == 1.] = targets_Z
	
	
	#==============================================================================
	# g.   D E L E T E   P A R T I C L E S   O U T S I D E   B O U N D A R I E S
	#==============================================================================
	# - deleting particles outside numerical boundaries 
	targets_R, targets_Z, active = delete_outer_particles(rminb, rmaxb, zminb, zmaxb, step_R, step_Z, targets_R, targets_Z, active)
	
	
	#==============================================================================
	# h.     D E L E T E   P A R T I C L E S   O U T S I D E   L I M I T E R
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
	# i.   C O M P U T E   D I S T A N C E S   B E T W E E N   P H I   S T E P S
	#==============================================================================
	# - connection length conditional
	if (lc_mode == 1 and tar_mode == 'line'):
	
		# - differentials and mean of radius between phi = k and phi = k-1
		delta_R =  targets_R_lc1 - targets_R_lc0
		delta_Z =  targets_Z_lc1 - targets_Z_lc0
		mean_R  = (targets_R_lc1 + targets_R_lc0) / 2 
		
		# - metric in cyllindrical coordinates between k and k-1 with element-wise
		# operations
		metric = np.sqrt( delta_R ** 2 + delta_Z ** 2 + mean_R ** 2 * delta_phi ** 2 )			
	
	
	#==============================================================================
	# j.           S U M M I N G   C O N N E C T I O N   L E N G T H
	#==============================================================================
	# - connection length conditional
	if (lc_mode == 1 and tar_mode == 'line'):
	
		# - summing to prior connection length for active particles only
		connection_length = connection_length + metric * active
	
	
	#==============================================================================
	# k.                 E X I T   I F   N O   T A R G E T S
	#==============================================================================
	# - exit loop if there are no more targets
	if len(targets_R) == 0:
		print('\n')
		print('=====================================================')	
		print('        ATTENTION: YOU RUN OUT OF PARTICLES!         ')
		print('=====================================================')
		break
		
		
	#==============================================================================
	# l.            S A V I N G   I N T E R S E C T I O N S
	#==============================================================================
	# - if full turn completed
	if (fraction == 0 and step > 0):

		# - step in process bar
		bar.next()
		
		# - append intersection
		intersect_R, intersect_Z = save_intersection(targets_R, targets_Z, intersect_R, intersect_Z)


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

# - declare number of intersections
saved_data.createDimension('R', len(intersect_R))

# - declare variables
plot_R = saved_data.createVariable('plot_R',  'f4', ('R',))#,  zlib=True)
plot_Z = saved_data.createVariable('plot_Z',  'f4', ('R',))#,  zlib=True)

# - assign values to variables
plot_R[:] = intersect_R
plot_Z[:] = intersect_Z

# - close netCDF data to be saved
saved_data.close()


#==============================================================================
# 17.        S A V I N G   C O N N E C T I O N   L E N G T H   D A T A   
#==============================================================================
# - connection length conditional
if (lc_mode == 1 and tar_mode == 'line'):

	# - Create a netCDF dataset to save data
	saved_data = nc.Dataset(str_save + '_Lc.nc', 'w')

	# - declare dimensions
	saved_data.createDimension('R', ntargets)

	# - declare variables
	points_R     = saved_data.createVariable('points_R',  'f4', ('R',))#,  zlib=True)
	points_Z     = saved_data.createVariable('points_Z',  'f4', ('R',))#,  zlib=True)
	connection_L = saved_data.createVariable('connection_L',  'f4', ('R',))#,  zlib=True)

	# - assign values to netcdf variables
	points_R[:]     = plot_targets_R
	points_Z[:]     = plot_targets_Z
	connection_L[:] = connection_length

	# - close netCDF data to be saved
	saved_data.close()


#==============================================================================
# 18.               P L O T   2 :   I N T E R S E C T I O N S
#==============================================================================
# - total number of intersections
N_int = len(intersect_R)

# - plotting intersections
plt.figure(2)
plt.plot(intersect_R, intersect_Z, 'ro', markersize=0.5)
plt.plot(plot_targets_R, plot_targets_Z, 'go', markersize=1.5)
plt.title('Field line intersections at $\phi$ = 0 - Total: ' + str(N_int))
plt.xlabel('R (m)')
plt.ylabel('Z (m)')


#==============================================================================
# 19.                  P L O T   2 :   L I M I T E R
#==============================================================================
# - conditional limiter mode
if lim_mode == 1:

	# - consider outer points as 1
	limiter_R, limiter_Z = limiter_boundary(grid_R, grid_Z, limiter[0,:,:])
	
	# - plot limiter line
	plt.plot(limiter_R, limiter_Z, 'ko', markersize=0.5)

# - save figure and show
plt.savefig(str_save + '.png', dpi=1000)
plt.show()


#==============================================================================
# 20.           P L O T   3 :   C O N N E C T I O N   L E N G T H
#==============================================================================
# - connection length conditional
if (lc_mode == 1 and tar_mode == 'line'):

	# - plotting connection length in vertical logarithmic scale
	plt.figure(2)
	plt.plot(plot_targets_R, connection_length, 'r-')
	plt.yscale('log')
	plt.title('Radial connection length')
	plt.xlabel('R (m)')
	plt.ylabel('L_c (m)')
	plt.savefig(str_save + '_Lc.png', dpi=1000)
	plt.show()


#==============================================================================
# 21.                          G O O D B Y E
#==============================================================================
print('\n')
print('=====================================================')
print('            FIELD LINE TRACING FINISHED!             ')
print('=====================================================')
print('\n')


#==============================================================================
#==============================================================================

