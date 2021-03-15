#==============================================================================
# PROGRAM - filit_v2.0.py
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
#> 15 03 2021 - Carlos Romero
#> 
#> TO DO:
#> - Comment 'plt.show()' for programmed shutdown
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

from matplotlib import pyplot as plt
from filit_mod  import read_data_hint,     \
					   read_limiter,       \
					   limiter_boundary,   \
					   create_targets,     \
					   evolve_cylindrical, \
					   inside_boundaries,  \
					   inside_limiter,     \
					   metric_cylindrical

from progress.bar import Bar

from scipy.interpolate import interpn, interp2d


#==============================================================================
# 1.                             I N P U T S
#==============================================================================	
tar_mode     = 'line-horizontal' # Target's distribution: 'line-vertical/horizontal' - 'matrix'

ntargets1D   =  2                # Number of targets per line

nturns       =  50#0              # Number of turns to the reactor

step_phi_mm  =  20                # Toroidal step in mm, for major radius R0

filepath     = 'snapfile.80.nc'  # Full name of the input file (HINT output)

limpath      = 'limiter.nc'      # Full name of the limiter/divertor geometry file

lim_mode     =  1                # Simulation boundaries (0), real limiter (1)

renormalize_dbg  =  1            # Boolean for renormalization of Br and Bz

dbg = 0                          # Continuity counter debug


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
print('     ····       ···  ·········· ···     ····   v2.0  ')
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
# ?.                D E B U G :   R E N O R M A L I Z A T I O N  
#==============================================================================
# - optional dynamic trajectory plot for debug
if renormalize_dbg == 1:
		 
	# - radial magnetic field
	B_R = B_R / 5

	# - vertical magnetic field
	B_Z = B_Z / 5


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
# 6.                  G R I D S   A N D   S T E P S  (R, Z)
#==============================================================================	
# - defining original grids
grid_phi = np.linspace(0,     2*np.pi, n_grid_phi)
grid_Z   = np.linspace(zminb, zmaxb,   n_grid_Z)
grid_R   = np.linspace(rminb, rmaxb,   n_grid_R)

# - computing steps of regular grids
step_Z = (zmaxb - zminb ) / (n_grid_Z - 1)
step_R = (rmaxb - rminb ) / (n_grid_R - 1)


#==============================================================================
# 7.     I N I T I A L   T A R G E T S   A N D   I N T E R S E C T I O N
#==============================================================================
# - define a radial left and right margins to create targets far from
# boundaries
margin_left   = 90  * step_R
margin_right  = 130 * step_R

# - analagous for vertical margins
margin_bottom = 70  * step_Z
margin_top    = 70  * step_Z

# - creating initial targets from given input values
initial_targets_R, initial_targets_Z = create_targets(tar_mode, ntargets1D, rminb + margin_left, \
														    rmaxb - margin_right,                \
															zminb + margin_bottom,               \
															zmaxb - margin_top)

# - depending on the selected target mode we have ntargets1D or ntargets1D^2
# total targets
ntargets = len(initial_targets_R)

# - declare array to save intersections
print('\n')
print('NOTE: saving intersections on symmetric toroidal planes! ')
# NOTE: multiplying by mtor allows us to save intersection points at every symmetry toroidal angle
intersect_R = np.nan * np.empty((ntargets, nturns * mtor))
intersect_Z = np.nan * np.empty((ntargets, nturns * mtor))

# - assing initial targets to 0-th column of intersections
intersect_R[:, 0] = np.copy(initial_targets_R)
intersect_Z[:, 0] = np.copy(initial_targets_Z)


#==============================================================================
# 8.           F O R   3-D   I N T E R P O L A T I O N   G R I D
#==============================================================================
# - create array of grid points in space
grid_points = (grid_phi, grid_Z, grid_R)


#==============================================================================
# 9.               F O R   C O N N E C T I O N   L E N G T H
#==============================================================================
# - initial definition of connection length 
connection_length = np.zeros(np.shape(intersect_R[:,0]))


#==============================================================================
# 10.                  M A I N   L O O P :   P A R T I C L E S
#==============================================================================
# - process bar for loop
print('\n')
bar_particles = Bar('Processing particles', max=ntargets)


# - particle loop
for particle in range(ntargets):

	bar_particles.next()
	
	#==============================================================================
	# *         P R E D E F I N I T I O N S   P E R   P A R T I C L E
	#==============================================================================
	# - starting point per particle
	dummy_phi0 = 0.
	dummy_Z0 = np.copy(initial_targets_Z[particle])	
	dummy_R0 = np.copy(initial_targets_R[particle])

	# - create array
	dummy0 = np.array([dummy_phi0, dummy_Z0, dummy_R0])
	
	# - intersection plane counter
	plane = 1
	
	# - turn counter
	turn = 1
	
		
	#==============================================================================
	# **            S E C O N D A R Y   L O O P :   T U R N S
	#==============================================================================
	# - while number of turns is not achieved
	while turn <= nturns:
	
	
		#==============================================================================
		# A.        3-D   I N T E R P O L A T I O N   G R I D :   F I E L D S
		#==============================================================================
		# - interpolate magnetic fields in 3D (it could also be done with the manually
		# defined interpolation functions (2D and 3D) in 'field_line_mod'
		B_phi_inter = interpn(grid_points, B_phi, dummy0, method='linear')
		B_Z_inter   = interpn(grid_points, B_Z,   dummy0, method='linear')
		B_R_inter   = interpn(grid_points, B_R,   dummy0, method='linear')
		
		# - create array
		B_inter = (B_phi_inter, B_Z_inter, B_R_inter)


		#==============================================================================
		# B.                T O R O I D A L   E V O L U T I O N
		#==============================================================================
		# - compute next point in field line 
		dummy1 = evolve_cylindrical(dummy0, step_phi_mm, B_inter)
		
				
		#==============================================================================
		# ?.              D E B U G :   F I E L D   C O N T I N U I T Y 
		#==============================================================================
		if dbg > 0:
			value = max(abs(B_inter - B_inter_copy))
			#print('value = ' + str(value))
			
			if value > 0.1:
				print('-------------------------')
				print('B_inter')
				print(B_inter)
				print('B_inter_copy')
				print(B_inter_copy)
				print('angle = ' + str(dummy1[0]))
				print('-------------------------')
			
		dbg = dbg + 1
		B_inter_copy = np.copy(B_inter)
		
		
		#==============================================================================
		# C.   D E L E T E   P A R T I C L E S   O U T S I D E   B O U N D A R I E S
		#==============================================================================
		# - deleting particles outside numerical boundaries 
		inside_bool = inside_boundaries(rminb, rmaxb, zminb, zmaxb, dummy1)
		
		
		#==============================================================================
		# D.     D E L E T E   P A R T I C L E S   O U T S I D E   L I M I T E R
		#==============================================================================
		# - conditional limiter mode
		if (inside_bool == 1 and lim_mode == 1):
		

			#==============================================================================
			# *        3-D   I N T E R P O L A T I O N   G R I D :   L I M I T E R
			#==============================================================================
			# - interpolate limiter
			limiter_inter = interpn(grid_points, limiter, dummy1)#, method='linear', bounds_error=True, fill_value=0.)
			
			# - delete particles outside limiter/divertor geometry
			# (we want to pass a value, NOT an array: that's the reason for '[0]')
			inside_bool = inside_limiter(limiter_inter[0], dummy1)
		
		
		#==============================================================================
		# E.                 E X I T   I F   N O   T A R G E T S
		#==============================================================================
		# - if particle is outside boundaries or outside limiter (optional)
		if inside_bool == 0:
		
			# - exit while loop if particle is outside -> go to next particle
			break
			
		
		#==============================================================================
		# F.   C O M P U T E   D I S T A N C E S   B E T W E E N   P H I   S T E P S
		#==============================================================================
		# - differential element of length between original and updated position of 
		# particle following a field line
		dlength = metric_cylindrical(dummy0, dummy1)
		
		# - add it to the total
		connection_length[particle] = connection_length[particle] + dlength

		
		#==============================================================================
		# G.            S A V I N G   I N T E R S E C T I O N S
		#==============================================================================
		# - the array components are:
		#angle0 = dummy0[0]
		#angle1 = dummy1[0]
		
		# - angles at which particle's intersection with toroidal plane will be saved
		symmetry_angle = 2 * np.pi / mtor
		
		# - boolean for plane crossing: if updated particle has passed the intersection
		# saving plane, then cross = 1; else cross = 0 (meaning they are in the same
		# toroidal sector)
		cross = (dummy1[0] % symmetry_angle < dummy0[0] % symmetry_angle) * 1
		
		# - update position
		dummy0 = np.copy(dummy1)
		
		# - delete array
		#np.delete(dummy1)
		
		
		# - if trajectory has crossed an intersection-saving plane
		if cross == 1:
		
			# - save positions
			intersect_Z[particle, plane] = np.copy(dummy0[1])
			intersect_R[particle, plane] = np.copy(dummy0[2])
			
			# - update intersection plane counter
			plane = plane + 1
			
			# - when we are back to the origin after 'mtor' symmetry planes
			if plane % mtor == 0:
			
				# - update turn
				turn = turn + 1
				


#==============================================================================
# 11.                E N D   M A I N   L O O P :   T U R N S
#==============================================================================	
# - ending bar counter and printing
bar_particles.finish()

#==============================================================================
# 12.                 D A T A   S A V I N G   P A T H 
#==============================================================================
# - folder name
folder = 'filit_results'

# - full string for saving data
str_save = folder + '/' + tar_mode + '_ntargets-' + str(ntargets) + '_nturns-' + str(nturns)


#==============================================================================
# 13.                S A V I N G   I N T E R S E C T I O N S   
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
# 14.        S A V I N G   C O N N E C T I O N   L E N G T H   D A T A   
#==============================================================================
# - declare dimensions
saved_data.createDimension('R', ntargets)

# - declare variables
plot_connection_length = saved_data.createVariable('plot_connection_length',  'f4', ('Np',))#,  zlib=True)

# - assign values to netcdf variables
plot_connection_length[:] = connection_length


#==============================================================================
# 15.         S A V I N G   L I M I T E R   B O U N D A R I E S
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
# 16.               P L O T   1 :   I N T E R S E C T I O N S
#==============================================================================
# - plotting intersections
plt.figure(2)

# - printing loop
for target in range(ntargets):

	plt.plot(intersect_R[target,:], intersect_Z[target,:], 'o', markersize=0.5)
	
plt.plot(initial_targets_R, initial_targets_Z, 'g*', markersize=1.5)
plt.title('Field line intersections at $\phi$ = 0')
plt.xlabel('R (m)')
plt.ylabel('Z (m)')
plt.grid(True)


#==============================================================================
# 17.                  P L O T   1 :   L I M I T E R
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
# 18.           P L O T   2 :   C O N N E C T I O N   L E N G T H
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
# 19.                          G O O D B Y E
#==============================================================================
print('\n')
print('=====================================================')
print('            FIELD LINE TRACING FINISHED!             ')
print('=====================================================')
print('\n')


#==============================================================================
#==============================================================================
