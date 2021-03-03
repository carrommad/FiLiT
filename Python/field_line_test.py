#------------------------------------------------------------------------------
# PROGRAM - field_line_test.py
#------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------
# Modules
import netCDF4           as nc
import numpy             as np
import timeit

from matplotlib     import pyplot as plt
from field_line_mod import save_intersection
						   
						   #------------------------------------------------------------------------------
#------------------------------------------------------------------------------			   
# EVERYTHING HAS BEEN SUCCESSFULLY TESTED IN:
#
# - field_line_mod
#
#------------------------------------------------------------------------------#------------------------------------------------------------------------------
