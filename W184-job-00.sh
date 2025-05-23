#!/bin/bash -l

# this job script is specifical for steps 02-04 of the pipeline that does not use
# more than two cores:
# one for the main thread and one worker to do the talys calculation
# so I will run it on a single node requsting two cores

#SBATCH -A naiss2025-22-537
#SBATCH -p main -n 1
#SBATCH -t 5:00
#SBATCH -J job-00-ND-pipeline

BASE_DIR=${HOME}/Public/NEPU
PROJ_DIR=/cfs/klemming/projects/supr/naiss2023-22-58
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-rackham-with-stdout-redirect-new.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/NEPU-scripts/script-Cr
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/NEPU-scripts/config/config-W184.R

module load PDC
module load apptainer

apptainer exec --bind $BASE_DIR --bind $PROJ_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/01_prepare_experimental_data.R $CONFIG_FILE
