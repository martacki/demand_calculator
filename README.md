# Interannual Electricity Demand Calculator

Large parts of this code were originally developed by [Lieke van der
Most](https://github.com/L-vdM) (University of Groningen) in the *EU renewable
energy modelling framework*. The original version of the code can be found
[here](https://github.com/L-vdM/EU-renewable-energy-modelling-framework) and is
referenced below as [1]. [1] has been validated against historical
electricity demand data reported on the [ENTSO-E transparancy
platform](https://transparency.entsoe.eu/). If you use this code, please refer to her work as well.

We have made the following adjustments to the original version:

- generate hourly instead of daily electricity consumption profiles
- use `snakemake` for workflow management
- trim repository to demand-related code and data
- adjust code to accept cutouts from `atlite` for weather data

## Purpose

Variations in weather conditions affect electricity demand patterns. This
workflow generates country-level electricity consumption time series based on
weather data using analysis by [Lieke van der Most](https://github.com/L-vdM)
correlating historical electricity demand to temperature. This workflow first
calculates a daily electricity demand based on the regression model developed in
[1]. Subsequently, cumulative daily electricity demands are disaggregated using
a hourly profile sampled from a random historical day (that is the same weekday)
from the [Open Power System
Database](https://data.open-power-system-data.org/time_series/). The resulting
`output/demand_hourly.csv` file is compatible with the open-source electricity
system model [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur).

Holidays are treated like weekend days. Data on national holidays across Europe
are obtained using another repository by [Aleksander
Grochowicz](https://github.com/aleks-g) and others that similarly computes
artificial electricity demand time series:
[github.com/aleks-g/multidecade-data](https://github.com/aleks-g/multidecade-data/blob/v1.0/load%20data/create_artificial_demand.ipynb).
The holidays are stored at `input_files/noworkday.csv`.

## Installation and Usage

### Clone the Repository

Download the [demand_calculator](https://github.com/martacki/demand_calculator)
repository using `git`.

```bash
/some/other/path % cd /some/path/without/spaces
/some/path/without/spaces % git clone https://github.com/martacki/demand_calculator.git
```

### Install Dependencies with conda/mamba

Use [`conda`](https://docs.conda.io/en/latest/miniconda.html) or
[`mamba`](https://github.com/QuantStack/mamba) to install the required packages
listed in
[environment.yaml](https://github.com/martacki/demand_calculator/blob/master/environment.yaml).

The environment can be installed and activated using

```bash
.../demand_calculator % conda env create -f environment.yaml
.../demand_calculator % conda activate demand
```

### Retrieve Input Data

The only required additional input files are ERA5 cutouts which can be recycled
from the [PyPSA-Eur weather data deposit on
Zenodo](https://zenodo.org/record/6382570#.Yx4KN2xByV4). Place the file
`europe-2013-era5.nc` in the following location (and rename!):

```
./input_files/cutouts/europe-era5-2013.nc
```

Cutouts for other weather years than 2013 can be built using the `build_cutout`
rule from the [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur) repository.

### Run the Workflow

This repository uses `snakemake` for workflow management. To run the complete
workflow, execute:

```bash
.../demand_calculator % snakemake -jall all
```

After successfully running the workflow, the output files will be located in
`output/energy_demand` named `demand_hourly_{yr}.csv`.

The years to compute can be modified directly in the `Snakefile`.

## Example Plots of the Workflow Output

An exemplary plot of the result for Germany (DE) of an exemplary week of January, 2013
comparing the results of this workflow with the [Open Power System Data](https://data.open-power-system-data.org/time_series/) in hourly resolution.

![ts-DE](https://user-images.githubusercontent.com/53824825/188666599-bff05561-e601-40d0-9e90-51a6eb68455c.png)

An exemplary plot of the result for Spain (ES) of an exemplary week of January, 2013
comparing the results of this workflow with the [Open Power System Data](https://data.open-power-system-data.org/time_series/) in hourly resolution.

![ts-ES](https://user-images.githubusercontent.com/53824825/188666633-9844a3d8-fc60-4940-ad57-eb92b61dd6a6.png)
