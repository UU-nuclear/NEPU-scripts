#!/bin/bash -l

#SBATCH -A naiss2023-22-58

#SBATCH -p node -N 3
#SBATCH -t 03:00:00
#SBATCH -J job-04-ND-pipeline

BASE_DIR=/proj/naiss2023-22-58
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-no-stdout-redirect-new.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/script-Cr
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/config/config-Cr52-mod-cov-func-try2.R

module load openmpi/4.0.2
mpirun -np 60 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/07_tune_talyspars_mod_cov_func.R $CONFIG_FILE
mpirun -np 2 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/07_5_addGPobs.R $CONFIG_FILE
mpirun -np 60 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/10_tune_talyspars_with_defect_mod_cov_func.R $CONFIG_FILE
