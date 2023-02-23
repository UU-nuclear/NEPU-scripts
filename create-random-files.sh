#!/bin/bash -l

# this job script is specifical for step 05 of the pipeline

#SBATCH -A naiss2023-22-58
#SBATCH -p node -N 5
#SBATCH -t 02:00:00
#SBATCH -J ND-random-files

BASE_DIR=/proj/naiss2023-22-58
SIF_FILE=$BASE_DIR/ND-eval-pipeline/NDeval-pipeline-rackham-with-stdout-redirect-new.sif
SCRIPT_DIR=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/script-Cr
CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/config/config-Cr52-mod-cov-func.R
#CONFIG_FILE=$BASE_DIR/ND-eval-pipeline/eval-fe56-scripts/config/config-test.R

module load openmpi/4.0.2
#mpirun -np 2 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/11_calculate_posterior_approximation_with_defect_GLS.R $CONFIG_FILE
mpirun -np 100 apptainer exec --bind $BASE_DIR $SIF_FILE Rscript --vanilla $SCRIPT_DIR/12_create_randomfiles_mod_cov_func.R $CONFIG_FILE
