#!/bin/bash -l

# this job script is specifical for step 02 of the pipeline that does not use
# more than two cores:
# one for the main thread and one worker to do the talys calculation
# so I will run it on a single node requsting two cores

#SBATCH -A naiss2024-22-324
#SBATCH -p core -n 2
#SBATCH -t 5:00
#SBATCH -J talysTemp

BASE_DIR=/proj/naiss2024-22-324
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-rackham-with-stdout-redirect.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/script-Cr
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/config/config-Cr52-mod-cov-func.R
#CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/config/config-test.R

module load openmpi/4.0.2
mpirun -np 2 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/07_tune_talyspars_mod_cov_func.R $CONFIG_FILE

