# SPDX-FileCopyrightText: : 2022 Martha Frysztacki (KIT), Fabian Neumann (TUB)
#
# SPDX-License-Identifier: MIT


rule all:
    input:
        expand(
            "output/energy_demand/demand_hourly_{yr}.csv", yr=list(range(1951, 2022))
        ),
    output:
        "output/energy_demand/demand_hourly.csv",
    resources:
        mem_mb=8000,
    run:
        import pandas as pd

        df = pd.concat([pd.read_csv(f, index_col=0) for f in input], axis=0)
        df.to_csv(output[0])


rule prepare_cutout:
    input:
        "input_files/cutouts/europe-era5-{yr}.nc",
    output:
        "climate_data/temperature_{yr}.nc",
    resources:
        mem_mb=2000,
    script:
        "model/prepare_cutout.py"


rule generate_daily_demand:
    input:
        population="input_files/demand_fit/population_t2m_grid.nc",
        shapefile_countries="input_files/EEZ_land_union_v3_202003/EEZ_Land_v3_202030.shp",
        demand_fit="input_files/demand_fit/demand_fit_values.nc",
        country_convertor="input_files/demand_fit/dict_population_per_country.json",
        noworkday="input_files/noworkday.csv",
        climate_data="climate_data/temperature_{yr}.nc",
    output:
        demand_daily="output/energy_demand/demand_daily_{yr}.nc",
        pop_per_country_nc="output/resources/pop_per_country_{yr}.nc",
        pop_per_country_json="output/resources/dict_pop_per_country_{yr}.nc",
        population_weights="output/resources/population_t2m_grid_weights_{yr}.nc",
        population="output/resources/population_t2m_grid_{yr}.nc",
    threads: 1
    resources:
        mem_mb=8000,
    script:
        "model/demand_daily.py"


rule generate_hourly_demand:
    input:
        demand_daily="output/energy_demand/demand_daily_{yr}.nc",
        reference=expand(
            "input_files/reference_demand/load_{yrs}.csv",
            yrs=[2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018],
        ),
    output:
        "output/energy_demand/demand_hourly_{yr}.csv",
    resources:
        mem_mb=2000,
    script:
        "model/demand_hourly.py"
