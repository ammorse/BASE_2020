#!/bin/sh
#SBATCH --mail-user=sherbina@usc.edu
#SBATCH --job-name=stan2_and_stan4c
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/mcintyre/share/ase_multi_env/scripts/slurm_logs/%x_%A_%a.out
#SBATCH -t 24:00:00
#SBATCH --ntasks=1              # single cpu core
#SBATCH --mem=3gb               # per processor memory

DIR_PROJECT=/blue/mcintyre/share/ase_multi_env

cd $DIR_PROJECT/programs

date;hostname;pwd

#Set the number of runs (i.e. number of columns in the design file ) that a
#SLURM task should run
PER_TASK=${1}

#Recalculate the starting and ending values for this array taks based on the
#SLURM task and the number of runs per task
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

echo -e "This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM"

module load R/3.6

DATA_FILE=${2}
qsim1=$(echo $DATA_FILE | rev | cut -d '/' -f1 | rev | cut -d'_' -f6)
qsim2=$(echo $DATA_FILE | rev | cut -d '/' -f1 | rev | cut -d'_' -f8)
MODEL_DESIGN_FILE=$DIR_PROJECT/design_files/${3}

for (( RUN=$START_NUM; RUN<=$END_NUM; RUN++ )); do
    echo $RUN
    line=$(tail -n +2 $MODEL_DESIGN_FILE | head -n $RUN | tail -n 1)
    STAN_PROGRAM=$(echo $line | awk -F',' '{print $1}')
    echo $STAN_PROGRAM
    ITER=$(echo $line | awk -F',' '{print $2}')
    BURNIN=$(echo $line  | awk -F',' '{print $3}')
    THIN=$(echo $line  | awk -F',' '{print $4}')
    REFRESH=$(echo $line  | awk -F',' '{print $5}')
    STEPSIZE=$(echo $line  | awk -F',' '{print $6}')
    are_q_equal=$(echo $line  | awk -F',' '{print $7}')
    if [[ $are_q_equal == Y ]]; then
      qmodel1=${qsim1}
      qmodel2=${qsim2}
    else
      echo '$are_q_equal is different from 1'
      qmodel1=$(echo $line  | awk -F',' '{print $8}')
      qmodel2=${qsim2}
    fi
    echo "are_q_equal" $are_q_equal 
    echo "qsim1" $qsim1 
    echo "qmodel1" $qmodel1 
    NCOND=$(echo $line  | awk -F',' '{print $9}')
    R --no-save --args data=${DATA_FILE} qmodel1=${qmodel1} qmodel2=${qmodel2} stan=${STAN_PROGRAM} iter=${ITER} burnin=${BURNIN} \
    thin=${THIN} refresh=${REFRESH} \
    stepsize=${STEPSIZE} ncond=${NCOND} < ${STAN_PROGRAM}_model_data_analysis_03kos.r
done
