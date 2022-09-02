# Interannual demand calculator (PyPSA-Eur compatible)

Large parts of this code were originally developed by our colleague Lieke van der Most (University of Groningen) in her EU renewable energy modelling framework. The original version of the code can be found here: https://github.com/L-vdM/EU-renewable-energy-modelling-framework.

The major differences are:
- generation of hourly instead of daily electricity consumption profiles for the given set of countries in Europe
- workflow handling - we use snakemake, such that the workflow is comparable to PyPSA-Eur

# Purpose

Variations in weather conditions typically imply variations in electricity demand patterns that are highly correlated.
This workflow provides a solution to generate electricity consumption time-series that depend on a given weather year and are compatible with the open source energy system model PyPSA-Eur (https://github.com/PyPSA/pypsa-eur).
