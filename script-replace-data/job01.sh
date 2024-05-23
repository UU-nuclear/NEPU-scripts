#!/bin/bash -l

#SBATCH -A naiss2024-22-324
#SBATCH -p core -n 10
#SBATCH -t 10:00
#SBATCH -J job-01-ND-pipeline

BASE_DIR=/proj/naiss2024-22-324
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-rackham-with-stdout-redirect-new.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/script-replace-data
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/config/test-replace-data.R

#module load openmpi/4.0.2
#mpirun -np 2 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/02_create_reference_calculation.R $CONFIG_FILE

#apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/03_extract_experimental_uncertainties.R $CONFIG_FILE

apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/04_replace_experimental_data.R $CONFIG_FILE

