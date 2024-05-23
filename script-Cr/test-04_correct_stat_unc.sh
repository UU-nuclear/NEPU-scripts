#!/bin/bash -l

#SBATCH -A naiss2024-22-324
#SBATCH -p core -n 20
#SBATCH -t 01:00:00
#SBATCH -J test04-ND-pipeline

BASE_DIR=/proj/naiss2024-22-324
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-rackham-with-stdout-redirect-new.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/script-Cr
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/config/config-Cr52-mod-cov-func-try4.R

apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/04_correct_stat_unc.R $CONFIG_FILE

