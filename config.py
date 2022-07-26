
import numpy as np
import glob

# timestep
dt = 24
# set namee of author and project for output files 
author_name = 'Martha Frysztacki'
project_name = 'calculate_demand'

# set regions
selected_countries = 'europe'
## set_area in lon, lat 
l = -31.25
r = 74.5
t = 79
b = 33

country_name_column ='TERRITORY1'
abbreviation_column = 'ISO_TER1'

### climate
runs = ['2013']
t2m_varname = 'temperature'
ro_varname = 'ro' # runoff
wind_varname = 'wind10m' # surface wind
t2mmax_varname = 't2mmax' # max daily temperature
rad_varname = 'rsds' # radiation 
dis_varname = 'discharge' # discharge