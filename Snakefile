# SPDX-FileCopyrightText: : 2022 Martha Frysztacki (KIT)
#
# SPDX-License-Identifier: MIT

rule prepare_cutout:
    input: "input_files/cutouts/europe-{yr}-era5.nc"
    output: "climate_data/temperature_{yr}.nc"
    resources: mem_mb=500
    script: "model/prepare_cutout.py"

rule generate_daily_demand:
    input:
        population          = "input_files/demand_fit/population_t2m_grid.nc",
        shapefile_countries = "input_files/EEZ_land_union_v3_202003/EEZ_Land_v3_202030.shp",
        demand_fit          = "input_files/demand_fit/demand_fit_values.nc",
        country_convertor   = "input_files/demand_fit/dict_population_per_country.json",
        climate_data        = "climate_data/temperature_{yr}.nc",
    output:
        demand_daily = "output/energy_demand/demand_daily_{yr}.nc"
    threads: 1
    resources: mem_mb=500
    script: "run.py"

rule generate_timeseries:
    input:
        demand_daily    = "output/energy_demand/demand_daily_{yr}.nc",
        reference       = expand("input_files/reference_demand/load_{yrs}.csv",
						         yrs=[2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018])
    output: "output/energy_demand/demand_hourly_{yr}.nc"
    script: "model/demand_timeseries.py"
