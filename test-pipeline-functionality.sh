#!/usr/bin/env bash

# BASE_DIR is the path that includes the script pipeline
BASE_DIR=/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/

# The path to the Apptainer image
SIF_FILE=$BASE_DIR/NEPU-with-rstudio.sif

# The path to where the EXFOR database is stored
EXFOR_DB_DIR=$BASE_DIR/../exfor_db/

# These paths are relative to BASE_DIR
SCRIPT_DIR=$BASE_DIR/eval-fe56-scripts/script-Cr
CONFIG_FILE=$BASE_DIR/eval-fe56-scripts/config/config-test-functionality.R


# start an instance of the pipeline apptainer/singularity image
# if the execution failed and the script did not finish you may need to stop
# this instance 'manually' by running the following command
# singularity instance stop NEPU
singularity instance start --bind $BASE_DIR --bind $EXFOR_DB_DIR $SIF_FILE NEPU

# start the exfor mongoDB
singularity exec instance://NEPU start_EXFOR_mongoDB.sh $EXFOR_DB_DIR

# retrieve the experimental data from exfor database, need to have the mongoDB running
singularity exec instance://NEPU Rscript --vanilla $SCRIPT_DIR/01_prepare_experimental_data.R $CONFIG_FILE

# do a first default talys calculation
mpirun -np 2 singularity exec instance://NEPU Rscript --vanilla $SCRIPT_DIR/02_create_reference_calculation.R $CONFIG_FILE

# run the extraction of the experimental uncertainties and the rule based correction of these
singularity exec instance://NEPU Rscript --vanilla $SCRIPT_DIR/03_extract_experimental_uncertainties.R $CONFIG_FILE

# 04_correct_stat_unc is skipped in the functionality test because it is very computationally heavy
# singularity exec instance://NEPU Rscript --vanilla $SCRIPT_DIR/04_correct_stat_unc.R $CONFIG_FILE

# run the Marginal Likelihood Optimization of the experimental uncertainties
singularity exec instance://NEPU Rscript --vanilla $SCRIPT_DIR/04_tune_experimental_uncertainties.R $CONFIG_FILE

# # calculate the full Jacobian matrix, this may take some time
mpirun -np 4 singularity exec  instance://NEPU Rscript --vanilla $SCRIPT_DIR/05_create_reference_jacobian.R $CONFIG_FILE

# using the full jacobian select only parameters affected by the experimental data for optimization
singularity exec instance://NEPU Rscript --vanilla $SCRIPT_DIR/05_select_parameters.R $CONFIG_FILE

# tune the hyper-parameters of the Gaussian Process on energy dependence of talys-parameters
singularity exec instance://NEPU Rscript --vanilla $SCRIPT_DIR/06_tune_endep_hyperpars_mod_cov_func.R $CONFIG_FILE

# perform the first LM parameter optimization (this will take time)
mpirun -np 4 singularity exec instance://NEPU  Rscript --vanilla $SCRIPT_DIR/07_tune_talyspars_mod_cov_func.R $CONFIG_FILE

# add the Gaussian Process in the observable
singularity exec instance://NEPU  Rscript --vanilla $SCRIPT_DIR/07_5_addGPobs.R $CONFIG_FILE

# perform the second LM parameter optimization (this will take time)
mpirun -np 4 singularity exec instance://NEPU  Rscript --vanilla $SCRIPT_DIR/10_tune_talyspars_with_defect_mod_cov_func.R $CONFIG_FILE

# calculate the posterior approximation
singularity exec instance://NEPU Rscript --vanilla $SCRIPT_DIR/11_calculate_posterior_approximation_with_defect_GLS.R $CONFIG_FILE

# generate random files from the posterior approximation
mpirun -np 4 singularity exec instance://NEPU Rscript --vanilla $SCRIPT_DIR/12_create_randomfiles_mod_cov_func.R $CONFIG_FILE

singularity instance stop NEPU