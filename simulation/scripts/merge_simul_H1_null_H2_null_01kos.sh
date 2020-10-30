#!/bin/bash

PROJ_DIR=/mnt/c/Users/ksherbina/mcintyre_lab/BASE_2020/simulation/programs
OUTPUT_DIR=${PROJ_DIR}/g3_sim_output/H1_null_H2_null_H3_null
python3 merge_2_simul_01kos.py H1_null 1 H1_null 2 ${OUTPUT_DIR}
