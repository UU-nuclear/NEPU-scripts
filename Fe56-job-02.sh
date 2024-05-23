#!/bin/bash -l

# this job script is specifical for step 05 of the pipeline

#SBATCH -A naiss2024-22-324

#SBATCH -p node -N 5
#SBATCH -t 01:20:00
#SBATCH -J job-02-ND-pipeline

BASE_DIR=/proj/naiss2024-22-324
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-rackham-with-stdout-redirect-new.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/script-Cr
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/config/config-Fe56.R

module load openmpi/4.0.2
mpirun -np 100 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/05_create_reference_jacobian.R $CONFIG_FILE
