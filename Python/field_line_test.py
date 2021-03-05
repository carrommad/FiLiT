#==============================================================================
# PROGRAM - field_line_test.py
#==============================================================================
#> 
#> DESCRIPTION:
#> Loads saved field line traces from 'field_line_tracer.py' in netCDF format.
#> 
#> REVISION HISTORY:
#> 26 02 2021 - Carlos Romero
#> 
#> TO DO:
#> -
#> 
#==============================================================================
#==============================================================================
# Modules
import netCDF4           as nc
import numpy             as np
import timeit

from matplotlib     import pyplot as plt
						   
#==============================================================================
#==============================================================================			   
# EVERYTHING HAS BEEN SUCCESSFULLY TESTED IN:
#
# - field_line_load   -> SUCCESSFULLY TESTED
# - field_line_mod    -> SUCCESSFULLY TESTED
# - field_line_tracer -> TESTING:



from scipy.interpolate import interpn


[[ 3.12932081  0.02847917  3.31939744]
 [ 3.12932081 ,-0.02717821  ,3.95221959]
 [ 3.12932081 ,-2.60669801 , 4.1232395 ]]




#
#==============================================================================
#==============================================================================
