# SPDX-FileCopyrightText: : 2022 Martha Frysztacki (KIT)
#
# SPDX-License-Identifier: MIT


rule generate_peakdemand:
    input:
        population          = "input_files/demand_fit/population_t2m_grid.nc",
        shapefile_countries = "input_files/EEZ_land_union_v3_202003/EEZ_Land_v3_202030.shp",
        demand_fit          = "input_files/demand_fit/demand_fit_values.nc",
        climate_data        = "input_files/climate_data/t2m_d/temperature_{yr}.nc",
        country_convertor   = "input_files/demand_fit/dict_population_per_country.json"
    output:
#        pop_per_country = "output/run3/energy_demand/pop_per_country_{yr}.nc",
        demand          = "output/run3/energy_demand/demand_{yr}.nc"
    threads: 1
    resources: mem_mb=500
    script: "run.py"