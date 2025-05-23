#!/bin/bash -l

# this job script is specific for steps 02-04 of the pipeline that does not use
# more than two cores:
# one for the main thread and one worker to do the talys calculation
# so I will run it on a single node requsting two cores

#SBATCH -A naiss2024-22-324
#SBATCH -p main -n 1
#SBATCH -t 5:00
#SBATCH -J job-05-ND-pipeline

BASE_DIR=${HOME}/Public/NEPU
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-rackham-with-stdout-redirect-new.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/NEPU-scripts/script-Cr
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/NEPU-scripts/config/config-Fe56.R

module load PDC
module load apptainer

apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/11_calculate_posterior_approximation_with_defect_GLS.R $CONFIG_FILE

