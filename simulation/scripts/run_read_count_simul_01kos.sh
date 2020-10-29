#!/bin/bash

### There are two arguments to sbatch aside from sbatch name:
### (1) Name of design file
### (2) Which set is it of simulations with this design file

sbatch simulate_read_counts_NBmodel_01kos.sbatch design_H1_null_01kos.csv 1
sbatch simulate_read_counts_NBmodel_01kos.sbatch design_H1_null_01kos.csv 2
sbatch simulate_read_counts_NBmodel_01kos.sbatch design_H1_not-null_01kos.csv 1
