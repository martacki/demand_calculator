# -*- coding: utf-8 -*-
# adapted from https://github.com/L-vdM/EU-renewable-energy-modelling-framework
import json
import sys

import numpy as np
import pandas as pd
import xarray as xr

sys.path.append("model/")

import constants.mappings as mappings
import data_processing.attributes as attributes
import data_processing.masking as masking
import energy_computation.demand as demand
from cdo import Cdo

cdo = Cdo()

if "snakemake" not in globals():
    from helper import mock_snakemake

    snakemake = mock_snakemake("generate_daily_demand", yr=2013)

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
print("COMPUTING DEMAND")

### make population weights
masking.cut_netcdf_into_regions(
    df_countries.EU_map,
    global_population,
    population_per_country_nc,
    population_per_country_json,
    shapefile_countries,
    country_indexes=df_countries.index_nr,
)
# regrid to temperature grid
cdo.remapsum(
    snakemake.input.climate_data,
    input=population_per_country_nc,
    output=population_tempgrid,
    readCdf=True,
    options="-f nc",
)

# take weights
pop_temp = xr.open_dataset(
    population_tempgrid,
)
weights = pop_temp / pop_temp.sum(dim=["lat", "lon"])
weights.to_netcdf(snakemake.output.population_weights)

### compute demand
fv = xr.open_dataset(fitvalues_file)
# loop over temperature files

r = snakemake.wildcards.yr
print("Computing demand for", r)

# get file names of run
climate_data = snakemake.input.climate_data
# open temperature data
ds_t2m = xr.open_dataset(climate_data)
# prepare temperature dataset
ds_t2m["temperature"] = ds_t2m["temperature"] - 273.15
ds_t2m["temperature"].attrs["units"] = "degC"

ds_t2m["temperature"].attrs.update(standard_name="temperature")
# prepare weighted temperature
ds_demand = xr.Dataset()
ds_demand["temp"] = (ds_t2m["temperature"] * weights.population).sum(
    dim=["lat", "lon"], keep_attrs=True
)
# to match dimensions of fitvalues (country,period) perfv weekend and weekday

country_mapping = (
    df_countries[["nuts_id", "index_nr"]].dropna().set_index("nuts_id").squeeze()
)
noworkday = pd.read_csv(snakemake.input.noworkday, index_col=0, parse_dates=True)
noworkday = noworkday.rename(columns=country_mapping.astype(float))
noworkday.columns.name = "country"
noworkday.index.name = "time"

df = ds_demand["temp"].to_pandas()

ds_demand = xr.concat(
    [
        xr.Dataset(dict(temp=df.where(~noworkday).stack().to_xarray())),
        xr.Dataset(dict(temp=df.where(noworkday).stack().to_xarray())),
    ],
    "period",
)
ds_demand["period"] = ["weekday", "weekend"]
# select only the countries that are in the input data and for which we have fitted data:
countries = list(
    set(np.unique(noworkday.columns)).intersection(set(np.unique(fv.country.data)))
)
ds_demand = ds_demand.sel(country=countries)
fv = fv.sel(country=countries)
# compute demand with fit variables
ds_demand = demand.compute_demand(ds_demand, fv)

ds_demand = ds_demand.sum(dim="period").transpose("time", "country")
ds_demand = attributes.set_global_attributes(
    ds_demand, "Entsoe-ERA5 fit and HW3", grid="gaussian n80", area="Europe"
)

# rename countries
with open(snakemake.input.country_convertor) as f:
    data = json.load(f)
ds_demand["country"] = [data[str(i)] for i in list(ds_demand.country.values)]

# save file
ds_demand.to_netcdf(snakemake.output[0])
