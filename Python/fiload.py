#==============================================================================
# PROGRAM - fiload.py
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
# Modules
import netCDF4  as nc
import numpy    as np

from matplotlib import pyplot as plt

#==============================================================================
# File path
filename = 'test'
filepath = filename + '.nc'

#==============================================================================
# Read data
dataset = nc.Dataset(filepath)

R  = dataset.variables['R_grid'][:]
Z  = dataset.variables['Z_grid'][:]
DD = dataset.variables['Density_matrix'][:]

#==============================================================================
# Plot 
plt.figure(1)
plt.pcolormesh(R, Z, DD, vmax=50)
plt.colorbar()
plt.title('Field line intersection density at $\phi$ = 0')
plt.xlabel('R (m)')
plt.ylabel('Z (m)')
plt.savefig(filename + '.png', dpi=1000)
plt.show()

#==============================================================================
#==============================================================================
