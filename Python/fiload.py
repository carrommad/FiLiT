#==============================================================================
# PROGRAM - field_line_load.py
#==============================================================================
#> 
#> DESCRIPTION:
#> Loads saved field line traces from 'field_line_tracer.py' in netCDF format.
#> 
#> REVISION HISTORY:
#> 23 02 2021 - Carlos Romero
#> 
#> TO DO:
#> -
#> 
#==============================================================================


#==============================================================================
# 0.                            M O D U L E S
#==============================================================================

import netCDF4  as nc
import numpy    as np

from matplotlib import pyplot as plt


#==============================================================================
# 1.                         F I L E   P A T H
#==============================================================================

filename = 'test'
filepath = filename + '.nc'

#==============================================================================
# 2.                         R E A D   D A T A
#==============================================================================

dataset = nc.Dataset(filepath)

plot_intersect_R  = dataset.variables['plot_intersect_R'][:]
plot_intersect_Z  = dataset.variables['plot_intersect_Z'][:]

plot_connection_length  = dataset.variables['plot_connection_length'][:]

plot_limiter_R  = dataset.variables['plot_limiter_R'][:]
plot_limiter_Z  = dataset.variables['plot_limiter_Z'][:]

ntargets = np.size(plot_intersect_R[:,0])

# - a way of seing if data is in a line or a matrix: basically this rules out
# vertical and horizontal lines
if plot_intersect_Z[0,0] != plot_intersect_Z[-1,0] and plot_intersect_R[0,0] != plot_intersect_R[-1,0]:
	ntargets1D = int(np.sqrt(ntargets))

#==============================================================================
# 4.               P L O T   1 :   I N T E R S E C T I O N S
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
#==============================================================================
