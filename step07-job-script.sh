#!/bin/bash -l

# this job script is specifical for step 05 of the pipeline

#SBATCH -A naiss2024-22-324

# I will try to ask for 2 nodes to start with
#SBATCH -p node -N 5
#SBATCH -t 02:00:00
#SBATCH -J ND-pipeline-step05

BASE_DIR=/proj/naiss2024-22-324
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-rackham-with-stdout-redirect.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/script-Cr
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/config/config-Cr52-mod-cov-func.R
#CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/config/config-test.R

module load openmpi/4.0.2
#mpirun -np 100 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/07_tune_talyspars_mod_cov_func.R $CONFIG_FILE
mpirun -np 2 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/07_5_addGPobs.R $CONFIG_FILE
mpirun -np 100 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/10_tune_talyspars_with_defect_mod_cov_func.R $CONFIG_FILE
