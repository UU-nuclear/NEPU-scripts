#!/bin/sh
# paramters:
#        $1 = the pipeline singularity .sif file
#        $2 = the location of the exfor DB
#        $3 = the pipeline configuration file (here config_mpi_reduced-pipeline.R)

singularity instance start $1 eval-inst
#singularity exec instance://eval-inst start_EXFOR_mongoDB.sh $2
#singularity exec instance://eval-inst Rscript --vanilla script/01_prepare_experimental_data.R $3
#singularity exec instance://eval-inst mpirun -np 1 Rscript --vanilla script/02_create_reference_calculation_mpi.R $3
#singularity exec instance://eval-inst Rscript --vanilla script/03_extract_reduced_experimental_uncertainties.R $3
#singularity exec instance://eval-inst Rscript --vanilla script/04_tune_experimental_uncertainties.R $3
#singularity exec instance://eval-inst mpirun -np 1 Rscript --vanilla script/05_create_reference_jacobian_mpi.R $3
#singularity exec instance://eval-inst Rscript --vanilla script/06_tune_endep_hyperpars.R $3

singularity exec instance://eval-inst mpirun -np 1 Rscript --vanilla script/07_tune_talyspars.R $3
singularity exec instance://eval-inst mpirun -np 1 Rscript --vanilla script/08_calculate_posterior_approximation.R $3
singularity exec instance://eval-inst mpirun -np 1 Rscript --vanilla script/09_create_randomfiles_mpi.R $3
singularity instance stop eval-inst

#run_step("01_prepare_experimental_data.R")
#run_step("02_create_reference_calculation.R")
#run_step("03_extract_reduced_experimental_uncertainties.R")
#run_step("04_tune_experimental_uncertainties.R")
#run_step("05_create_reference_jacobian.R")
#run_step("06_tune_endep_hyperpars.R")
#run_step("07_tune_talyspars.R")
#run_step("08_calculate_posterior_approximation.R")
#run_step("09_create_randomfiles.R")
