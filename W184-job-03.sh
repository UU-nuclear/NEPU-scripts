#!/bin/bash -l

# this job script is the pipeline script that does not use
# more than a single node:

#SBATCH -A naiss2025-22-537
#SBATCH -p main -n 20
#SBATCH -t 03:00:00
#SBATCH -J job-03-ND-pipeline

BASE_DIR=${HOME}/Public/NEPU
PROJ_DIR=/cfs/klemming/projects/supr/naiss2023-22-58
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-rackham-with-stdout-redirect-new.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/NEPU-scripts/script-Cr
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/NEPU-scripts/config/config-W184.R

module load PDC
module load apptainer

apptainer exec --bind $BASE_DIR --bind $PROJ_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/05_select_parameters.R $CONFIG_FILE
apptainer exec --bind $BASE_DIR --bind $PROJ_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/06_tune_endep_hyperpars_mod_cov_func.R $CONFIG_FILE

# for the case in config-Cr52-mod-cov-func.R the optimization took about 1 h on 20 cores, so 20 corehours

