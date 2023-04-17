#!/bin/bash -l

# this job script is specifical for steps 02-04 of the pipeline that does not use
# more than two cores:
# one for the main thread and one worker to do the talys calculation
# so I will run it on a single node requsting two cores

#SBATCH -A naiss2023-22-58
#SBATCH -p core -n 20
#SBATCH -t 1:00
#SBATCH -J job-01-ND-pipeline

BASE_DIR=/proj/naiss2023-22-58
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-rackham-with-stdout-redirect-new.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/script-Cr
SCRIPT_DIR_ALT=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/script
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/config/config-Fe56.R

apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR_ALT/03_extract_experimental_uncertainties.R $CONFIG_FILE

apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/04_correct_stat_unc.R $CONFIG_FILE

apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/04_tune_experimental_uncertainties.R $CONFIG_FILE

