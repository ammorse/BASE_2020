#!/bin/bash
sbatch --array=1-2 simulate_and_stan_model_data_analysis_02fm.sbatch 2 design_H1_not_null_H2_not_null_H3_null_01kos.csv design_stan_model_g1_qmodel_and_qsim_same_01kos.csv
sbatch --array=1-2 simulate_and_stan_model_data_analysis_02fm.sbatch 2 design_H1_not_null_01fm.csv design_stan_model_1env_01fm.csv
