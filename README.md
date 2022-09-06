# Interannual demand calculator (PyPSA-Eur compatible)

Large parts of this code were originally developed by our colleague Lieke van der Most (University of Groningen) in her EU renewable energy modelling framework. The original version of the code can be found here: https://github.com/L-vdM/EU-renewable-energy-modelling-framework and is referenced below as [1].

The major differences are:
- generation of hourly instead of daily electricity consumption profiles for the given set of countries in Europe
- workflow handling - we use snakemake, such that the workflow is comparable to PyPSA-Eur

# Purpose

Variations in weather conditions typically are highly correlated with variations in electricity demand patterns and thus imply variations.
This workflow provides a solution to generate electricity consumption time-series that depend on a given weather year and are compatible with the open source energy system model [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur).
This workflow calculates a daily electricity demand based on an evaluated package [1]. Then, it dis-aggregates the cumulative daily electricity demand to an hourly profile by sampling a random historical day (that is the same week-day) from the Open Power System Database (https://data.open-power-system-data.org/time_series/). The daily electricity demand is then distributed proportional to this sampled histrical data, resulting in an hourly profile. The resulting .csv document is compatible with PyPSA-Eur.

# How To ...

We set up the workflow to be easy to handle for [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur) users by using the same workflow managment system `snakemake`. The only required input file, an `era5` cutout, can be recycled from [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur) as well. The workflow automatically converts it to the required format.

The [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur) cutout as has to be placed in ``input_files/cutouts``. After running the snakemake workflow, the output ``.csv`` file will be located in ``output/energy_demand/``.

## Installation

### Clone the Repository

First of all, clone the [demand_calculator](https://github.com/martacki/demand_calculator) repository using the version control system ``git``.
The path to the directory into which the ``git repository`` is cloned, must **not** have any spaces! If you do not have ``git`` installed, follow installation instructions [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).

    /some/other/path % cd /some/path/without/spaces

    /some/path/without/spaces % git clone https://github.com/martacki/demand_calculator.git

### Install Python Dependencies

The workflow relies on a set of other Python packages to function.
We recommend using the package manager and environment management system ``conda`` to install them.
Install [miniconda](https://docs.conda.io/en/latest/miniconda.html), which is a mini version of [Anaconda](https://www.anaconda.com/) that includes only ``conda`` and its dependencies or make sure ``conda`` is already installed on your system.
For instructions for your operating system follow the ``conda`` [installation guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

The python package requirements are curated in the [envs/environment.yaml](https://github.com/martacki/demand_calculator/blob/master/environment.yaml) file. The environment can be installed and activated using

    .../pypsa-eur % conda env create -f envs/environment.yaml

    .../pypsa-eur % conda activate pypsa-eur

Note that activation is local to the currently open shell!
After opening a new terminal window, one needs to reissue the second command!

If you have troubles with a slow ``conda`` installation, we recommend to install [mamba](https://github.com/QuantStack/mamba) as a fast drop-in replacement via

    conda install -c conda-forge mamba

and then install the environment with

    mamba env create -f envs/environment.yaml

