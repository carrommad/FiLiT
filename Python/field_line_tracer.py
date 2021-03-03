#------------------------------------------------------------------------------
# PROGRAM - field_line_tracer.py
#------------------------------------------------------------------------------
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
#> 03 03 2021 - Carlos Romero
#> 
#> TO DO:
#> - What if R and Z dimensions are different?
#> -
#> 
#------------------------------------------------------------------------------
# Modules
import netCDF4           as nc
import numpy             as np

import timeit
import time

from matplotlib     import pyplot as plt
from field_line_mod import read_data_hint, create_targets, grid_neighbours,     \
						   interpolate_field, step_phi, delete_outer_particles, \
						   save_intersection

from progress.bar import Bar


#------------------------------------------------------------------------------
# INPUTS
#------------------------------------------------------------------------------		
 
mode       = 'matrix'          # Target's distribution: 'line' - 'matrix'

ntargets1D =  3                # Number of targets per line

nturns     =  100              # Number of turns to the reactor

filepath   = 'snapfile.80.nc'  # Full name of the input file (HINT output)

debug      =  0                # Debug mode: inactive(0), active (1)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Welcome and terminal output
print('\n')
print('===============================')
print(' Welcome to Field Line Tracer! ')
print('===============================')

print('\n')
print('=====================================================')
print(' You have selected:                                  ')
print('   - full input file name       -> ' + filepath       )
print('   - target mode: line/matrix   -> ' + mode           )
print('   - number of targets per line -> ' + str(ntargets1D))
print('   - number of turns            -> ' + str(nturns)    )
print('=====================================================')
	
# - debug conditional printing
if debug == 1:
	pause =  0.1

	print('\n')
	print('===================================')
	print(' ATTENTION: running on debug mode! ')
	print('===================================')
	print(' You have selected pause = ' + str(pause) + ' s')

#------------------------------------------------------------------------------
# Read data
dataset = nc.Dataset(filepath)
mtor, rminb, rmaxb, zminb, zmaxb, B_R, B_phi, B_Z = read_data_hint(dataset)

print('\n')
print('=========================')
print(' Data successfully read! ')
print('=========================')
print('\n')

# - debug conditional printing
if debug == 1:
	print('\n')
	print('======================================')
	print(' DEBUG: printing from read_data_hint! ')
	print(' - [rminb, rmaxb] =')
	print('[' + str(rminb) +', ' + str(rmaxb) + ']')
	time.sleep(pause)
	print(' - [zminb, zmaxb] =')
	print('[' + str(zminb) +', ' + str(zmaxb) + ']')
	time.sleep(pause)
	print(' - Dimensions B_R =')
	print(B_R.shape)
	time.sleep(pause)
	print(' - Dimensions B_phi =')
	print(B_phi.shape)
	time.sleep(pause)
	print(' - Dimensions B_Z =')
	print(B_Z.shape)
	time.sleep(pause)
	print('======================================')

#------------------------------------------------------------------------------
# Extracting dimensions from data and creating original grid
grid_size  = B_R.shape

# - careful with the order!!! I am not even sure about this one!!!
# -> IMPORTANT FOR interpolate_field
n_grid_phi = grid_size[0]
n_grid_Z   = grid_size[1]
n_grid_R   = grid_size[2]

grid_R = np.linspace(rminb,rmaxb,n_grid_R)
grid_Z = np.linspace(zminb,zmaxb,n_grid_R)

step_R = (rmaxb - rminb ) / (n_grid_R - 1)
step_Z = (zmaxb - zminb ) / (n_grid_R - 1)

# - debug conditional printing
if debug == 1:
	print('\n')
	print('==============================================')
	print(' DEBUG: printing recreation of original grid! ')
	print(' - Dimensions grid_R:')
	print(grid_R.shape)
	time.sleep(pause)
	print(' - Dimensions grid_Z:')
	print(grid_Z.shape)
	time.sleep(pause)
	print(' - [step_R, step_Z] =')
	print('[' + str(step_R) +', ' + str(step_Z) + ']')
	time.sleep(pause)
	print('==============================================')
	
#------------------------------------------------------------------------------
# Creating initital target points
# - These are read as inputs: mode = 'matrix' - ntargets1d = 100
targets_R, targets_Z = create_targets(mode, ntargets1D, rminb + 2*step_R, \
														rmaxb - 2*step_R, \
														zminb + 2*step_Z, \
														zmaxb - 2*step_Z)
													  
# - Depending on the selected mode we have ntargets1D or ntargets1D^2 total targets
ntargets = len(targets_R)

# - In case you are an idiot in the previous step: delete particles outside boundaries
targets_R, targets_Z = delete_outer_particles(rminb, rmaxb, zminb, zmaxb, step_R, step_Z, targets_R, targets_Z)

# - debug conditional printing
if debug == 1:
	print('\n')
	print('==========================================')
	print(' DEBUG: printing targets after filtering! ')
	print(' - Number of 1D targets    = ' + str(ntargets1D))
	time.sleep(pause)
	print(' - Number of total targets = ' + str(ntargets))
	time.sleep(pause)
	print(' - Dimensions targets_R    = ')
	print(targets_R.shape)
	time.sleep(pause)
	print(' - Dimensions targets_Z    = ')
	print(targets_Z.shape)
	time.sleep(pause)
	print('==========================================')

#------------------------------------------------------------------------------
# Iterate phi a certain number of toroidal turns and save at every turn
# - This is read as input: nturns = 1000
nsteps = nturns * n_grid_phi
intersect_R = targets_R
intersect_Z = targets_Z

# - debug conditional counting
if debug == 1:
	bar = Bar('Processing steps', max=nsteps)
	
	print('\n')
	print('=======================================')
	print(' DEBUG: printing info inside loop! ')
	print(' - Number of total targets = ' + str(len(targets_R)))

elif debug == 0:
	bar = Bar('Processing turns', max=nturns)


#------------------------------------------------------------------------------
# MAIN LOOP
#------------------------------------------------------------------------------
for step in range(nsteps):
		
	# - Fraction of the angle 
	angle = step % n_grid_phi
	
	# - Compute the grid neighbours of our targets
	index_R, index_Z, dist_R, dist_Z = grid_neighbours(targets_R, targets_Z, grid_R, grid_Z)

	# + Radial
	B_R_inter   = interpolate_field(index_R, index_Z, dist_R, dist_Z, B_R[angle,:,:], step_R, step_Z)
	# + Toroidal
	B_phi_inter = interpolate_field(index_R, index_Z, dist_R, dist_Z, B_phi[angle,:,:], step_R, step_Z)
	# + Vertical
	B_Z_inter   = interpolate_field(index_R, index_Z, dist_R, dist_Z, B_Z[angle,:,:], step_R, step_Z)
	
	# - Evolve finite differences
	delta_phi = 2 * np.pi / (n_grid_phi - 1) / mtor
	targets_R, targets_Z = step_phi(B_R_inter, B_phi_inter, B_Z_inter, delta_phi, rminb, rmaxb, targets_R, targets_Z)
	
	# - Delete particles outside boundaries
	targets_R, targets_Z = delete_outer_particles(rminb, rmaxb, zminb, zmaxb, step_R, step_Z, targets_R, targets_Z)
	
	# - Exit loop if there are no more targets
	if len(targets_R) == 0:
		print('\n')
		print('======================================')
		print(' ATTENTION: you run out of particles! ')
		print('======================================')
		break
		
	# - debug conditional counting steps
	if debug == 1:
		bar.next()
			
		#print(' - Number of total targets at step ' + str(step+1) + '/' + str(nsteps) + ' = ' + str(len(targets_R)))
		print(' - Number of total targets = ' + str(len(targets_R)))
		time.sleep(pause)
		
	# - Save intersection data at every turn of the stellarator/tokamak
	if angle == 0:
		
		# - debug conditional counting turns
		if debug == 0:
			bar.next()
		
		# - append intersection
		intersect_R, intersect_Z = save_intersection(targets_R, targets_Z, intersect_R, intersect_Z)
#------------------------------------------------------------------------------
# END OF MAIN LOOP
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# - ending bar counter and printing
bar.finish()

if debug == 1:
	print('==============================================')

#------------------------------------------------------------------------------
# Create a new thinner grid
nr_plot = 1 * n_grid_R
nz_plot = 1 * n_grid_Z

grid_R_update = np.linspace(rminb, rmaxb, nr_plot)
grid_Z_update = np.linspace(zminb, zmaxb, nz_plot)

density = np.zeros((nr_plot,nz_plot))

#------------------------------------------------------------------------------
# Compute and save intersection density from all intersected particles
for m in range(len(intersect_R)):	

	i = np.argmin(abs(grid_R_update - intersect_R[m]))
	j = np.argmin(abs(grid_Z_update - intersect_Z[m]))
	
	density[i,j] = density[i,j] + 1

#------------------------------------------------------------------------------
# Declare name for file to be saved
str_save = mode + '_ntargets-' + str(ntargets) + '_nturns-' + str(nturns)
	
#------------------------------------------------------------------------------
# Create a netCDF dataset to save data
saved_data = nc.Dataset(str_save + '.nc', 'w')

# - declare dimensions
saved_data.createDimension('R', nr_plot)
saved_data.createDimension('Z', nz_plot)

# - declare variables
R_grid = saved_data.createVariable('R_grid',  'f4', ('R',))#,  zlib=True)
Z_grid = saved_data.createVariable('Z_grid',  'f4', ('Z',))#,  zlib=True)
Density_matrix = saved_data.createVariable('Density_matrix', 'f4', ('R','Z',))#, zlib=True)

# - assign values to variables
R_grid[:] = grid_R_update
Z_grid[:] = grid_Z_update
Density_matrix[:] = density

print('\n')
print('=============================')
print(' Intersection density saved! ')
print('=============================')

#------------------------------------------------------------------------------
# Plot intersections
N_int = len(intersect_R)

# - impact plot deactivated
if 0:
	plt.figure(2)
	plt.plot(intersect_R, intersect_Z,'.r')
	plt.savefig('points.png')
	plt.show()

plt.figure(1)
plt.pcolormesh(grid_R_update, grid_Z_update, density)#, vmax=50)
plt.colorbar()
plt.title('Field line intersection density at $\phi$ = 0 - Total: ' + str(N_int))
plt.xlabel('R (m)')
plt.ylabel('Z (m)')
plt.savefig(str_save + '.png')
plt.show()

#------------------------------------------------------------------------------
# Closing netCDF data to be saved
saved_data.close()

# Print goodbye
print('\n')
print('=============================')
print(' Field Line Tracer finished! ')
print('=============================')
print('\n')

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
