#!/bin/bash -l

# this job script is specifical for step 05 of the pipeline

#SBATCH -A naiss2024-22-324

#SBATCH -p main -N 1
#SBATCH -t 01:20:00
#SBATCH -J job-02-ND-pipeline

BASE_DIR=${HOME}/Public/NEPU
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-rackham-with-stdout-redirect-new.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/NEPU-scripts/script-Cr
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/NEPU-scripts/config/config-Fe56.R

module load PDC
module load apptainer
module load openmpi/4.1.2-gcc12.2.0

mpirun -np 128 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/05_create_reference_jacobian.R $CONFIG_FILE
