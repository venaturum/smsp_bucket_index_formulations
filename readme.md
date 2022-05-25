# Bucket Index Formulations for the Single Machine Scheduling Problem

The code in this repository provides implementations (with either PuLP or gurobipy) for Bucket Index (BI) formulations, and Time Index formulation for the Single Machine Scheduling Problem.

The BI formulations are taken from 

[Clement, Riley. "MIXED INTEGER LINEAR PROGRAMMING MODELS FOR MACHINE SCHEDULING." (2015).](https://hdl.handle.net/1959.13/1310105)


### Install

Available as a python package `smsp_bi` which can be installed with pip:

    pip install git+https://github.com/venaturum/smsp_bucket_index_formulations

or with poetry (place the following dependency in pyproject.toml)

    smsp_bi = { git = "https://github.com/venaturum/smsp_bucket_index_formulations.git" }

Example usage can be found in *example.ipynb*.

### TODO

- Add cut separation algorithm
- Add BI-3 formulation
- Add BI-n formulation
