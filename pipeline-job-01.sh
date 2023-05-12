#!/bin/bash -l

# this job script is specifical for steps 02-04 of the pipeline that does not use
# more than two cores:
# one for the main thread and one worker to do the talys calculation
# so I will run it on a single node requsting two cores

#SBATCH -A naiss2023-22-58
#SBATCH -p core -n 2
#SBATCH -t 5:00
#SBATCH -J job-01-ND-pipeline

BASE_DIR=/proj/naiss2023-22-58
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-rackham-with-stdout-redirect-new.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/script-Cr
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/config/config-Cr52-hetGP.R

module load openmpi/4.0.2
mpirun -np 2 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/02_create_reference_calculation.R $CONFIG_FILE

apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/03_extract_experimental_uncertainties.R $CONFIG_FILE


