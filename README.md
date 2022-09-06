# Interannual demand calculator (PyPSA-Eur compatible)

Large parts of this code were originally developed by our colleague Lieke van der Most (University of Groningen) in her EU renewable energy modelling framework. The original version of the code can be found here: https://github.com/L-vdM/EU-renewable-energy-modelling-framework and is referenced below as [1].

The major differences are:
- generation of hourly instead of daily electricity consumption profiles for the given set of countries in Europe
- workflow handling - we use snakemake, such that the workflow is comparable to PyPSA-Eur

# Purpose

Variations in weather conditions typically are highly correlated with variations in electricity demand patterns and thus imply variations.
This workflow provides a solution to generate electricity consumption time-series that depend on a given weather year and are compatible with the open source energy system model PyPSA-Eur (https://github.com/PyPSA/pypsa-eur).
This workflow calculates a daily electricity demand based on an evaluated package [1]. Then, it dis-aggregates the cumulative daily electricity demand to an hourly profile by sampling a random historical day (that is the same week-day) from the Open Power System Database (https://data.open-power-system-data.org/time_series/). The daily electricity demand is then distributed proportional to this sampled histrical data, resulting in an hourly profile. The resulting .csv document is compatible with PyPSA-Eur.

# How To ...

We set up the workflow to be easy to handle for `PyPSA-Eur` users by using the same workflow managment system `snakemake`. The only required input file, an `era5` cutout, can be recycled from `PyPSA-Eur` as well. The workflow automatically converts it to the required format.

The `PyPSA-Eur` cutout as has to be placed in `input_files/cutouts`. After running the snakemake workflow, the output `.csv` file will be located in `output/energy_demand/`.
