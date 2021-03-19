#==============================================================================
# PROGRAM - field_line_load.py
#==============================================================================
#> 
#> DESCRIPTION:
#> Loads saved field line traces from 'field_line_tracer.py' in netCDF format.
#> 
#> REVISION HISTORY:
#> 17 03 2021 - Carlos Romero
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

from filib      import plot_Poincare_Lc


#==============================================================================
# 1.                           I N P U T S
#==============================================================================
# - name of the file
filename = 'results/matrix_ntargets-10000_nturns-100'

# - extension of the file
filepath = filename + '.nc'

# - boolean for plotting initial target distribution
plot_initial = 0

#==============================================================================
# 2.                         R E A D   D A T A
#==============================================================================
# - import netcdf data
dataset = nc.Dataset(filepath)

# - coordinates of intersections with \phi=0
intersect_R = dataset.variables['plot_intersect_R'][:]
intersect_Z = dataset.variables['plot_intersect_Z'][:]

# - connection length
connection_length = dataset.variables['plot_connection_length'][:]

# - limiter data may be optional
limiter_Z  = dataset.variables['plot_limiter_Z'][:]
limiter_R  = dataset.variables['plot_limiter_R'][:]


#==============================================================================
# 3.  P L O T :   P O I N C A R É   A N D   C O N N E C T I O N   L E N G T H
#==============================================================================
# - call predefined function
plot_Poincare_Lc(filename, intersect_Z, intersect_R, connection_length, limiter_Z, limiter_R, plot_initial)


#==============================================================================
#==============================================================================

