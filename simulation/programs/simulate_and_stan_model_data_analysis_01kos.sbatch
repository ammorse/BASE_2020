#!/bin/sh
#SBATCH --mail-user=fabio.marroni@uniud.it
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH --ntasks=1              # single cpu core
#SBATCH --mem=3gb               # per processor memory

DIR_PROJECT=/blue/mcintyre/fabio.marroni/2020_simul_stan

cd $DIR_PROJECT/programs

date;hostname;pwd

#Set the number of runs (i.e. number of columns in the design file ) that a
#SLURM task should run
PER_TASK=${1}

#Recalculate the starting and ending values for this array taks based on the
#SLURM task and the number of runs per task
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

echo -e "This is task $SLURM_ARAY_TASK_ID, which will do runs $START_NUM to $END_NUM"

module load R/3.6

SIMUL_DESIGN_FILE=$DIR_PROJECT/design_files/${2}
MODEL_DESIGN_FILE=$DIR_PROJECT/design_files/${3}

for (( RUN=$START_NUM; RUN<=$END_NUM; RUN++ )); do
  SCENARIO=$(tail -n +2 $SIMUL_DESIGN_FILE | head -n $RUN | tail -n 1)
  READS=$(echo $SCENARIO | awk -F',' '{print $1}')
  MULT=$(echo $SCENARIO | awk -F',' '{print $2}')
  qsim1=$(echo $SCENARIO | awk -F',' '{print $3}')
  qsim2=$(echo $SCENARIO | awk -F',' '{print $4}')
  alpha=$(echo $SCENARIO | awk -F',' '{print $5}')
  delta=$(echo $SCENARIO | awk -F',' '{print $6}')
  simruns=$(echo $SCENARIO | awk -F',' '{print $7}')
  if [ ${alpha} -eq 1 ]; then H1="null"; else H1="not_null"; fi
  if [ ${delta} -eq 1 ]; then H2="null"; else H2="not_null"; fi
  if [ ${alpha} -eq ${delta} ]; then H3="null"; else H3="not_null"; fi
  OUT_DIR=$DIR_PROJECT/g3_sim_output/H1_${H1}_H2_${H2}_H3_${H3}
  mkdir -p $OUT_DIR
  OUT_FILE=out_alpha_${alpha}_delta_${delta}_qsim_g1_${qsim1}_qsim_g2_${qsim1}_mult_${MULT}_myreads_${READS}_simruns_${simruns}.csv
  R --no-save --args alpha=${alpha} delta=${delta} q1=${qsim1} q2=${qsim2} mult=${MULT} reads=${READS} simruns=${SIMRUNS} outfile=${OUT_DIR}/${OUT_FILE} < simulate_read_counts_NBmodel_01kos.r
  tail -n +2 $MODEL_DESIGN_FILE | while IFS= read -r line; do
    STAN_PROGRAM=$(echo $line | awk -F',' '{print $1}')
    ITER=$(echo $line | awk -F',' '{print $2}')
    BURNIN=$(echo $line  | awk -F',' '{print $3}')
    THIN=$(echo $line  | awk -F',' '{print $4}')
    REFRESH=$(echo $line  | awk -F',' '{print $5}')
    STEPSIZE=$(echo $line  | awk -F',' '{print $6}')
    are_q_equal=$(echo $line  | awk -F',' '{print $7}')
    if [ $are_q_equal -eq 1 ]; then
      qmodel1=${qsim1}
      qmodel2=${qsim2}
    else
      qmodel1=$(echo $line  | awk -F',' '{print $8}')
      qmodel2=${qsim2}
    fi
    R --no-save --args data=${OUT_DIR}/${OUT_FILE} qmodel1=${qmodel1} qmodel2=${qmodel2} stan=${STAN_PROGRAM} iter=${ITER} burnin=${BURNIN} thin=${THIN} refresh=${REFRESH} \
    stepsize=${STEPSIZE} < ${STAN_PROGRAM}_model_data_analysis_02fm.r
  done
done
