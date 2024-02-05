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
#> RECOMMENDATIONS:
#> 1) Faster for high numbers of particles and low number of turns.
#>
#> - If magnetic periodicity ('Btor') is reduced from its original value 'mtor'
#> interpolation will be less accurate.
#> -
#> 
#> REVISION HISTORY:
#> 22 01 2024 - Carlos Romero
#> 
#> TO DO:
#> - Comment 'plt.show()' for programmed shutdown.
#> - In RK4, when doing interpolation: careful with 'fill_value=np.nan'!
#> - Write more comprehensively 'step_phi_cylindrical_rk4'!
#> - Try to write a fast version of the code without unnecessary copies and 
#> transition arrays!
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
				       step_phi_cylindrical_rk4,       \
				       delete_outer_particles,         \
				       delete_outer_particles_limiter, \
				       plot_Poincare_Lc

from progress.bar import Bar

from scipy.interpolate import interpn

#import IPython


#==============================================================================
# 1.                             I N P U T S
#==============================================================================	
tar_mode     = 'matrix'           # Target's distribution: 'line-horizontal/vertical' - 'matrix'

ntargets1D   =  100               # Number of targets per line (tip: use more targets and less turns)

nturns       =  1000              # Number of turns to the reactor

step_phi_mm  =  80                # Toroidal step in mm, for major radius R0

input_file   = 'snapfile.80.nc'   # Full name of the input file (HINT output)

lim_path     = 'limiter.nc'       # Full name of the limiter/divertor geometry file (it can be empty: '')


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
print('     ····       ···  ·········· ···     ····   v1.2  ')
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
print('   - full input file name       ->  ' + input_file       )
print('   - full limiter file name     ->  ' + lim_path       )
print('=====================================================')


#==============================================================================
# 3.                      R E A D   D A T A   H I N T
#==============================================================================	
# - magnetic equilibrium dataset
dataset = nc.Dataset(input_file)

# - read and extract data
mtor, rminb, rmaxb, zminb, zmaxb, B_phi_original, B_Z_original, B_R_original = read_data_hint(dataset)

# - make copy of original data for possible replication
B_phi = np.copy(B_phi_original)
B_Z   = np.copy(B_Z_original)
B_R   = np.copy(B_R_original)

# - define magnetic toroidal symmetry
Btor = 10
print('\n')
print('NOTE: selected magnetic periodicity: Btor = ' + str(Btor))

# - replicate symmetric toroidal data to have a copy of Btor-times the original,
# completing the stellarator toroidally
for n in range(Btor-1):
	B_phi = np.append(B_phi, B_phi_original, axis = 0)
	B_Z   = np.append(B_Z,   B_Z_original,   axis = 0)
	B_R   = np.append(B_R,   B_R_original,   axis = 0)
	
# - to solve periodicity problems: ngrid + 1
B_phi = np.append(B_phi, [B_phi_original[0,:,:]], axis = 0)
B_Z   = np.append(B_Z,   [B_Z_original[0,:,:]],   axis = 0)
B_R   = np.append(B_R,   [B_R_original[0,:,:]],   axis = 0)

# - close dataset
dataset.close()
	

#==============================================================================
# 4.                   R E A D   D A T A   L I M I T E R  
#==============================================================================	
# - define limiter mode depending on existence of limiter path
if len(lim_path) != 0:
	lim_mode = 1           
else:
	lim_mode = 0
	
# - limiter conditional: simulation boundaries (0), real limiter (1) 
if lim_mode == 1:

	# - limiter/divertor geometry data
	limset = nc.Dataset(lim_path)
	
	# - read and extract data
	limiter_original = read_limiter(limset)
	
	# - make copy of original data for possible replication
	limiter = np.copy(limiter_original)
	
	# - replicate symmetric toroidal data 'mtor - 1' number of times, to have
	# a copy of mtor-times the original, completing the stellarator toroidally
	for n in range(mtor-1):
		limiter = np.append(limiter, limiter_original, axis = 0)
		
	# - to solve periodicity problems: ngrid + 1
	limiter = np.append(limiter, [limiter_original[0,:,:]], axis = 0)
	
	# - close dataset
	limset.close()
	
	#==============================================================================
	# a.         L I M I T E R   D I M E N S I O N S   A N D   G R I D
	#==============================================================================	
	# - extracting dimensions from data
	grid_size_lim  = limiter.shape

	# - recreating original grid: careful with the order of original grids!
	# - original order [r,z,phi] but it appears transposed: [phi,z,r]
	n_grid_phi_lim = grid_size_lim[0] - 1

	
	#==============================================================================
	# b.          L I M I T E R   G R I D S   A N D   S T E P S  (R, Z)
	#==============================================================================	
	# - recreating original grid of periodic dimension
	grid_phi_lim = np.linspace(0, 2*np.pi, n_grid_phi_lim + 1)


#==============================================================================
# 5.         M A G N E T I C   D I M E N S I O N S   A N D   G R I D
#==============================================================================	
# - extracting dimensions from data
grid_size_B  = B_R.shape

# - recreating original grid: careful with the order
# 	+ original order [r,z,phi] 
#   + data appears transposed: [phi,z,r]
# - substraction of the added slice to periodic grid
n_grid_phi_B = grid_size_B[0] - 1
n_grid_Z     = grid_size_B[1]
n_grid_R     = grid_size_B[2]


#==============================================================================
# 6.   M A G N E T I C   G R I D   D E N S I T Y   F R O M   S T E P  (phi)
#==============================================================================	
# - computing major radius 
R0 = (rminb + rmaxb) / 2

# - for a given toroidal step, compute grid density only if bigger than
# original
n_grid_phi_B_enhanced = max(n_grid_phi_B, int(2 * np.pi * R0 * 1000 / step_phi_mm ))

# - ad hoc magnetic grid density
if True:
	print('\n')
	print('NOTE: ad hoc magnetic grid density! ')
	n_grid_phi_B_enhanced = n_grid_phi_lim

# - conditional for a too big toroidal step
if n_grid_phi_B_enhanced == n_grid_phi_B:
	print('\n')
	print('=====================================================')
	print('  ATTENTION: USING ORIGINAL MAGNETIC GRID SIZE: ' + str(n_grid_phi_B))
	print('=====================================================')

# - printing grid density info if step size actually enhances original grid
else:
	print('\n')
	print('=====================================================')
	print('    INFO: ENHANCED GRID vs. ORIGINAL: ' + str(n_grid_phi_B_enhanced) + ' > ' + str(n_grid_phi_B))
	print('=====================================================')


#==============================================================================
# 7.           M A G N E T I C   G R I D S   A N D   S T E P S  (R, Z)
#==============================================================================	
# - recreating original grid of periodic dimension
grid_phi_B = np.linspace(0, 2*np.pi, n_grid_phi_B + 1)

# - recreating original grids of non-periodic dimensions
grid_Z   = np.linspace(zminb, zmaxb,   n_grid_Z)
grid_R   = np.linspace(rminb, rmaxb,   n_grid_R)

# - we provide a different name to the toroidal step as it is also a
# differential for the connection length computation and to avoid confusion 
# with the imported function 'step_phi_cylindrical_rk4'
# - due to periodicity, the grid size is differently defined
delta_phi = 2 * np.pi / n_grid_phi_B_enhanced


# - computing steps of regular grids
step_Z = (zmaxb - zminb ) / (n_grid_Z - 1)
step_R = (rmaxb - rminb ) / (n_grid_R - 1)


#==============================================================================
# 8.                      I N I T I A L   T A R G E T S
#==============================================================================
# - define a radial left and right margins to create targets far from
# boundaries
left  = rminb + 30 * step_R
right = rmaxb - 30 * step_R

# - ad hoc target values
if False:
	print('\n')
	print('NOTE: ad hoc target horizontal positions! ')
	left  = 3.23
	right = 3.74

# - analagous for vertical margins
bottom = zminb + 30 * step_Z
top    = zmaxb - 30 * step_Z

# - creating targets for given input values
targets_Z, targets_R = create_targets(tar_mode, ntargets1D, left, right, bottom, top)

# - depending on the selected target mode we have ntargets1D or ntargets1D^2
# total targets
ntargets = len(targets_R)


#==============================================================================
# 9.   N U M B E R   O F   S T E P S   A N D   I N T E R S E C T I O N
#==============================================================================
# - iterate phi a certain number of steps given by grid and toroidal turn
nsteps = nturns * n_grid_phi_B_enhanced + 1

# - nplanes factor: number of planes per turn where intersection data will be
# saved
nplanes = 1 # mtor
print('\n')
print('NOTE: planes per turn where intersections are saved: nplanes = ' + str(nplanes))

# - number of planes to save intersection
nplanes_total  = nturns * nplanes + 1

# - declare array to save intersections
intersect_Z = np.nan * np.empty((ntargets, nplanes_total))
intersect_R = np.nan * np.empty((ntargets, nplanes_total))

# - assing initial targets to 0-th column of intersections
intersect_Z[:, 0] = np.copy(targets_Z)
intersect_R[:, 0] = np.copy(targets_R)


#==============================================================================
# 10.           F O R   3-D   I N T E R P O L A T I O N   G R I D
#==============================================================================
# - create array of grid points in space for magnetic field
grid_points_3D_B = (grid_phi_B, grid_Z, grid_R)

# - create array of grid points in space for limiter
grid_points_3D_lim = (grid_phi_lim, grid_Z, grid_R)


#==============================================================================
# 11.              F O R   C O N N E C T I O N   L E N G T H
#==============================================================================
# - initial definition of connection length 
connection_length = np.zeros(np.shape(targets_R))


#==============================================================================
# 12.              F O R   P A R T I C L E   T R A C K I N G
#==============================================================================
# - full copy of targets at prior iteration step in the 'for' loop
targets_Z_prior = np.copy(targets_Z)
targets_R_prior = np.copy(targets_R)

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
	targets_Z_prior[active == 1.] = np.copy(targets_Z)
	targets_R_prior[active == 1.] = np.copy(targets_R)

	
	#==============================================================================
	# b.             T O R O I D A L   S T E P   A N D   A N G L E
	#==============================================================================
	# - fraction of the angle 
	fraction = step % n_grid_phi_B_enhanced
	
	# - angle in radians modulo 2*pi
	angle = fraction * delta_phi
	
	
	#==============================================================================
	# c.             E V O L V E :   R U N G E - K U T T A    4  
	#==============================================================================
	# - create array of phi position of targets at this step
	targets_phi = angle * np.ones(np.shape(targets_Z))
	
	# - evolve finite differences with Runge-Kutta of 4th order
	targets_Z, targets_R = step_phi_cylindrical_rk4(grid_points_3D_B, B_phi, B_Z, B_R, targets_phi, targets_Z, targets_R, delta_phi)
	
	
	#==============================================================================
	# d.   D E L E T E   P A R T I C L E S   O U T S I D E   B O U N D A R I E S
	#==============================================================================
	# - deleting particles outside numerical boundaries 
	targets_Z, targets_R, active = delete_outer_particles(rminb, rmaxb, zminb, zmaxb, targets_Z, targets_R, active)
	
	
	#==============================================================================
	# e.     D E L E T E   P A R T I C L E S   O U T S I D E   L I M I T E R
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
		limiter_inter = interpn(grid_points_3D_lim, limiter, targets, method='linear', bounds_error=True, fill_value=0.)
		
		# - delete particles outside limiter/divertor geometry
		targets_Z, targets_R, active = delete_outer_particles_limiter(limiter_inter, targets_Z, targets_R, active)
		
		
	#==============================================================================
	# f.   C O M P U T E   D I S T A N C E S   B E T W E E N   P H I   S T E P S
	#==============================================================================
	# - differentials and mean of radius between phi = k and phi = k-1
	delta_Z =  targets_Z - targets_Z_prior[active == 1.]
	delta_R =  targets_R - targets_R_prior[active == 1.]
	mean_R  = (targets_R + targets_R_prior[active == 1.]) / 2 
	
	# - metric in cyllindrical coordinates between k and k-1 with element-wise
	# operations
	metric = np.sqrt( delta_R ** 2 + delta_Z ** 2 + mean_R ** 2 * delta_phi ** 2 )			

	# - summing to prior connection length for active particles only
	connection_length[active == 1.] = connection_length[active == 1.] + metric

	
	#==============================================================================
	# g.                 E X I T   I F   N O   T A R G E T S
	#==============================================================================
	# - if there are no more targets
	if len(targets_R) == 0:
		print('\n')
		print('=====================================================')	
		print('        ATTENTION: YOU RUN OUT OF PARTICLES!         ')
		print('=====================================================')
		
		# - exit loop 
		break
		
		
	#==============================================================================
	# h.            S A V I N G   I N T E R S E C T I O N S
	#==============================================================================
	# - saving at selected nplanes planes	
	if (fraction % (n_grid_phi_B_enhanced / nplanes) == 0 and step > 0):
		
		# - save positions
		intersect_Z[active == 1., plane] = np.copy(targets_Z)
		intersect_R[active == 1., plane] = np.copy(targets_R)

		# - update intersection plane counter
		plane = plane + 1
		
		# - if full turn completed
		if fraction == 0:
				
			# - step in process bar
			bar.next()
		

#==============================================================================
# 14.              E N D   M A I N   L O O P :   T U R N S
#==============================================================================	
# - ending bar counter and printing
bar.finish()


#==============================================================================
# 15.                 D A T A   S A V I N G   P A T H 
#==============================================================================
# - folder name
folder = ''
#folder = 'results/'

# - full string for saving data
output_file = folder + tar_mode + '_ntargets-' + str(ntargets) + '_nturns-' + str(nturns) + '_Btor-' + str(Btor)


#==============================================================================
# 16.                S A V I N G   I N T E R S E C T I O N S   
#==============================================================================
# - create a netCDF dataset to save data
saved_data = nc.Dataset(output_file + '.nc', 'w')

# - declare number of targets
saved_data.createDimension('Nt', ntargets)

saved_data.createDimension('Nt1D', ntargets1D)

# - declare number of planes where the intersections where saved
saved_data.createDimension('Np', nplanes_total)

# - declare variables
plot_intersect_Z = saved_data.createVariable('plot_intersect_Z',  'f4', ('Nt','Np',))#,  zlib=True)
plot_intersect_R = saved_data.createVariable('plot_intersect_R',  'f4', ('Nt','Np',))#,  zlib=True)

plot_targets_Z = saved_data.createVariable('plot_targets_Z',  'f4', ('Nt1D',))#,  zlib=True)
plot_targets_R = saved_data.createVariable('plot_targets_R',  'f4', ('Nt1D',))#,  zlib=True)

# - assign values to variables
plot_intersect_Z[:] = intersect_Z
plot_intersect_R[:] = intersect_R


#==============================================================================
# 17.        S A V I N G   C O N N E C T I O N   L E N G T H   D A T A   
#==============================================================================
# - declare dimensions
saved_data.createDimension('R', ntargets)

# - declare variables
plot_connection_length = saved_data.createVariable('plot_connection_length',  'f4', ('Nt',))#,  zlib=True)

# - assign values to netcdf variables
plot_connection_length[:] = connection_length


#==============================================================================
# 18.         S A V I N G   L I M I T E R   B O U N D A R I E S
#==============================================================================
# - conditional limiter mode
if lim_mode == 1:

	# - consider outer points as 1
	limiter_Z, limiter_R = limiter_boundary(grid_Z, grid_R, limiter[0,:,:])

# - if there is no limiter
elif lim_mode == 0:
	
	# - save an empty list
	limiter_Z = []
	limiter_R = []
	
# - declare number of intersections
saved_data.createDimension('L', len(limiter_R))
	
# - declare variables
plot_limiter_Z = saved_data.createVariable('plot_limiter_Z',  'f4', ('L',))#,  zlib=True)
plot_limiter_R = saved_data.createVariable('plot_limiter_R',  'f4', ('L',))#,  zlib=True)

# - assign values to variables
plot_limiter_Z[:] = limiter_Z
plot_limiter_R[:] = limiter_R

# - close netCDF data to be saved
saved_data.close()


#==============================================================================
# 19. P L O T :   P O I N C A R É   A N D   C O N N E C T I O N   L E N G T H
#==============================================================================
# - boolean for plotting initial target distribution and limiter
plot_initial = 0
plot_limiter = 0

# - call predefined function
plot_Poincare_Lc(output_file, targets_Z, intersect_Z, intersect_R, connection_length, limiter_Z, limiter_R, plot_initial, plot_limiter)


#==============================================================================
# 20.                          G O O D B Y E
#==============================================================================
print('\n')
print('=====================================================')
print('            FIELD LINE TRACING FINISHED!             ')
print('=====================================================')
print('\n')


#==============================================================================
#==============================================================================

