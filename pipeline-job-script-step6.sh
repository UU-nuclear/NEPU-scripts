#!/bin/bash -l

# this job script is the pipeline script that does not use
# more than a single node:

#SBATCH -A naiss2023-22-58
#SBATCH -p core -n 8
#SBATCH -t 00:15:00
#SBATCH -J talysTemp

BASE_DIR=/proj/naiss2023-22-58
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-no-stdout-redirect.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/script-Cr
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/config/config-Cr52-mod-cov-func.R
#CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/config/config-test.R

#apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/05_select_parameters.R $CONFIG_FILE
apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/06_tune_endep_hyperpars_mod_cov_func.R $CONFIG_FILE

