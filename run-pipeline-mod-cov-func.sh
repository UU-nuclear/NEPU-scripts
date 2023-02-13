#!/bin/sh
# paramters:
#        $1 = the pipeline singularity .sif file
#        $2 = the location of the exfor DB
#        $3 = the pipeline configuration file
# 
apptainer instance start --bind /TMC/alf $1 pipeline-inst # -- bind /TMC/alf will mount this directory in the container

#echo "**********************************************************************"
#echo "*                        step01                                      *"
#echo "**********************************************************************"
#apptainer exec instance://pipeline-inst start_EXFOR_mongoDB.sh $2
#apptainer exec instance://pipeline-inst Rscript --vanilla script-Cr/01_prepare_experimental_data.R $3

echo "**********************************************************************"
echo "*                        step02                                      *"
echo "**********************************************************************"
apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script-Cr/02_create_reference_calculation.R $3

#echo "**********************************************************************"
#echo "*                        step03                                      *"
#echo "**********************************************************************"
#apptainer exec instance://pipeline-inst Rscript --vanilla script-Cr/03_extract_experimental_uncertainties.R $3
#apptainer exec instance://pipeline-inst Rscript --vanilla script/visualization/plotExpData_withUnc.R $3
#
#echo "**********************************************************************"
#echo "*                        step04                                      *"
#echo "**********************************************************************"
#apptainer exec instance://pipeline-inst Rscript --vanilla script-Cr/04_tune_experimental_uncertainties.R $3
#
#echo "**********************************************************************"
#echo "*                        step05                                      *"
#echo "**********************************************************************"
#apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script-Cr/05_create_reference_jacobian.R $3
#
#echo "**********************************************************************"
#echo "*                        step06                                      *"
#echo "**********************************************************************"
#apptainer exec instance://pipeline-inst Rscript --vanilla script-Cr/06_tune_endep_hyperpars_mod_cov_func.R $3
#
#echo "**********************************************************************"
#echo "*                        step07                                      *"
#echo "**********************************************************************"
#apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script-Cr/07_tune_talyspars_mod_cov_func.R $3 
#
#echo "**********************************************************************"
#echo "*                        step07.5                                    *"
#echo "**********************************************************************"
#apptainer exec instance://pipeline-inst Rscript --vanilla script-Cr/07_5_addGPobs.R $3
#
#echo "**********************************************************************"
#echo "*                        step10                                      *"
#echo "**********************************************************************"
#apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script-Cr/10_tune_talyspars_with_defect_mod_cov_func.R $3
#
#echo "**********************************************************************"
#echo "*                        step11                                      *"
#echo "**********************************************************************"
#apptainer exec instance://pipeline-inst Rscript --vanilla script-Cr/11_calculate_posterior_approximation_with_defect_GLS.R $3
#
#echo "**********************************************************************"
#echo "*                        step12                                      *"
#echo "**********************************************************************"
#apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script-Cr/12_create_randomfiles_mod_cov_func.R $3

apptainer instance stop pipeline-inst