#!/bin/sh
#SBATCH --mail-user=sherbina@usc.edu
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --job-name=simulate_read_counts
#SBATCH --output=/blue/mcintyre/share/ase_multi_env/scripts/slurm_logs/%x_%j.out
#SBATCH --ntasks=1              # single cpu core
#SBATCH --mem=2gb               # per processor memory
#SBATCH --time=1:00:00

DIR_PROJECT=/blue/mcintyre/share/ase_multi_env

cd $DIR_PROJECT/programs

SIMUL_DESIGN_FILE=$DIR_PROJECT/design_files/${1}

module load R/3.6

for line in $(tail -n +2 $SIMUL_DESIGN_FILE); do
  READS=$(echo $line | awk -F',' '{print $1}')
  MULT=$(echo $line | awk -F',' '{print $2}')
  qsim1=$(echo $line | awk -F',' '{print $3}')
  qsim2=$(echo $line | awk -F',' '{print $4}')
  alpha=$(echo $line | awk -F',' '{print $5}')
  SIMRUNS=$(echo $line | awk -F',' '{print $6}')
  NCOND=$(echo $line | awk -F',' '{print $7}')
  echo ${NCOND}
  if [ ${alpha} == 1 ]; then H1="null"; else H1="not_null"; fi
  OUT_DIR=$DIR_PROJECT/g3_sim_output/H1_${H1}
  echo $OUT_DIR
  OUT_FILE=out_set_${2}_alpha_${alpha}_qsim_g1_${qsim1}_qsim_g2_${qsim1}_mult_${MULT}_myreads_${READS}_simruns_${SIMRUNS}.csv
  echo $OUT_FILE
  R --no-save --args alpha=${alpha} q1=${qsim1} q2=${qsim2} mult=${MULT} reads=${READS} simruns=${SIMRUNS} outfile=${OUT_DIR}/${OUT_FILE} ncond=${NCOND} < simulate_read_counts_NBmodel_03kos.r
done
