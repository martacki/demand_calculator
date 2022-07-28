import sys
import xarray as xr
import numpy as np
import json

sys.path.append('model/')

import constants.mappings as mappings
import data_processing.masking as masking
import data_processing.attributes as attributes
import energy_computation.demand as demand

from cdo import Cdo
cdo = Cdo()

df_countries = mappings.df_countries_select
# demand uses EU_map so drop nan values 
df_countries = df_countries.loc[df_countries.EU_map.notna()]
global_population = snakemake.input.population
shapefile_countries = snakemake.input.shapefile_countries
fitvalues_file = snakemake.input.demand_fit
population_tempgrid = snakemake.output.population
population_per_country_nc = snakemake.output.pop_per_country_nc
population_per_country_json = snakemake.output.pop_per_country_json

# =============================================================================
# Compute demand
# =============================================================================
print('COMPUTING DEMAND')

### make population weights
masking.cut_netcdf_into_regions(df_countries.EU_map, global_population, 
                                population_per_country_nc, population_per_country_json, shapefile_countries, 
                                country_indexes=df_countries.index_nr )
# regrid to temperature grid
cdo.remapsum(snakemake.input.climate_data,
     input=population_per_country_nc,
     output=population_tempgrid,
     readCdf=True,
     options='-f nc')
# take weights 
pop_temp = xr.open_dataset(population_tempgrid,)
weights = pop_temp/pop_temp.sum(dim=['lat', 'lon'])
weights.to_netcdf(snakemake.output.population_weights)

### compute demand
fv = xr.open_dataset(fitvalues_file)
# loop over temperature files

r = snakemake.wildcards.yr
print('Computing demand for', r)

# get file names of run
climate_data = snakemake.input.climate_data
# open temperature data
ds_t2m = xr.open_dataset(climate_data)
# prepare temperature dataset
ds_t2m["temperature"] = ds_t2m["temperature"] - 273.15
ds_t2m["temperature"].attrs['units'] = 'degC'

ds_t2m["temperature"].attrs.update(standard_name = "temperature")
# prepare weighted temperature
ds_demand = xr.Dataset()
# next row kills the job
ds_demand['temp'] = (ds_t2m["temperature"] * weights.population).sum(
    dim=['lat','lon'], keep_attrs=True)
# to match dimensions of fitvalues (country,period) per weekend and weekday
ds_demand = xr.concat(
    [ds_demand.where(ds_demand['time.dayofweek']<5, drop=True),
     ds_demand.where(ds_demand['time.dayofweek']>=5, drop=True)],
    'period')
ds_demand['period'] = ['weekday', 'weekend']
# select only the countries that are in the input data and for which we have fitted data:
countries = list(set(np.unique(ds_demand.country)).intersection(set(np.unique(fv.country.data))))
ds_demand = ds_demand.sel(country = countries)
fv = fv.sel(country = countries)
# compute demand with fit variables
ds_demand = demand.compute_demand(ds_demand, fv)
# remove period dimension
ds_demand = xr.concat([ds_demand.isel(period=0), ds_demand.isel(period=1)], 'time')
# update attributes
ds_demand = attributes.set_global_attributes(
    ds_demand, 'Entsoe-ERA5 fit and HW3', grid='gaussian n80', area='Europe')
    
# clean  file
ds_demand['demand'] = ds_demand.demand.transpose('time', 'country')
ds_demand = ds_demand.drop('period').dropna(dim='time')
# rename countries
with open(snakemake.input.country_convertor) as f:
    data = json.load(f)
ds_demand['country'] = [data[str(i)] for i in list(ds_demand.country.values)]
    
# save file
ds_demand.to_netcdf(snakemake.output[0])
