#!/bin/bash -l

#SBATCH -A naiss2023-22-58

#SBATCH -p node -N 3
#SBATCH -t 01:00:00
#SBATCH -J job-04-global

BASE_DIR=/proj/naiss2023-22-58
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-rackham-with-stdout-redirect-new.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/script-global-fit
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/config/config-Cr52-global-fit.R

module load openmpi/4.0.2
mpirun -np 60 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/07_tune_talyspars_mod_cov_func.R $CONFIG_FILE
