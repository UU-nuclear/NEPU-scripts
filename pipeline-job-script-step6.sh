#!/bin/bash -l

# this job script is the pipeline script that does not use
# more than a single node:

#SBATCH -A naiss2024-22-324
#SBATCH -p core -n 20
#SBATCH -t 02:00:00
#SBATCH -J talysTemp

BASE_DIR=/proj/naiss2024-22-324
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-rackham-with-stdout-redirect.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/NEPU-scripts/script-Cr
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/NEPU-scripts/config/config-Cr52-mod-cov-func.R
#CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/NEPU-scripts/config/config-test.R

#apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/05_select_parameters.R $CONFIG_FILE
apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/06_tune_endep_hyperpars_mod_cov_func.R $CONFIG_FILE

# for the case in config-Cr52-mod-cov-func.R the optimization took about 1 h on 20 cores, so 20 corehours

