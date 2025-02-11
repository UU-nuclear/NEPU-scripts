#!/bin/bash -l

#SBATCH -A naiss2024-22-324

#SBATCH -p node -N 2
#SBATCH -t 01:00:00
#SBATCH -J job-04-global

BASE_DIR=/proj/naiss2024-22-324
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-rackham-with-stdout-redirect-new.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/NEPU-scripts/script-global-fit
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/NEPU-scripts/config/config-Cr52-global-fit.R

module load openmpi/4.0.2
mpirun -np 40 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/07_tune_talyspars_mod_cov_func.R $CONFIG_FILE
