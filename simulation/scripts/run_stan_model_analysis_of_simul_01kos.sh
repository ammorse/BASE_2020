#!/bin/bash

SIMUL=/blue/mcintyre/share/ase_multi_env/g3_sim_output/H1_null_H2_null_H3_null/alpha1_1_alpha2_1_qsim-g1_0.8_qsim-g2_0.8_mult_1_myreads_100_simruns_1000.csv
sbatch --array=1-2 stan_model_data_analysis_01kos.sbatch 1 ${SIMUL} design_stan_model_g1_qmodel_and_qsim_same_01kos.csv
