#!/bin/bash -l

#SBATCH -A naiss2024-22-324

#SBATCH -p node -N 3
#SBATCH -t 03:00:00
#SBATCH -J job-04-ND-pipeline

BASE_DIR=/proj/naiss2024-22-324
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-rackham-with-stdout-redirect-new.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/NEPU-scripts/script-Cr
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/NEPU-scripts/config/config-Fe56.R

module load openmpi/4.0.2
mpirun -np 60 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/07_tune_talyspars_mod_cov_func.R $CONFIG_FILE
mpirun -np 2 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/07_5_addGPobs.R $CONFIG_FILE
mpirun -np 60 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/10_tune_talyspars_with_defect_mod_cov_func.R $CONFIG_FILE
