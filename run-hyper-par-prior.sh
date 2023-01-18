#!/bin/sh
# paramters:
#        $1 = the pipeline singularity .sif file
#        $2 = the location of the exfor DB
#        $3 = the pipeline configuration file

apptainer instance start --bind /TMC/alf $1 pipeline-inst # -- bind /TMC/alf will mount this directory in the container
# apptainer exec instance://pipeline-inst start_EXFOR_mongoDB.sh $2
# echo "**********************************************************************"
# echo "*                        step01                                      *"
# echo "**********************************************************************"
# apptainer exec instance://pipeline-inst Rscript --vanilla script/01_prepare_experimental_data.R $3
# echo "**********************************************************************"
# echo "*                        step02                                      *"
# echo "**********************************************************************"
# apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script/02_create_reference_calculation_with_unc_extra_endep_pars.R $3
# echo "**********************************************************************"
# echo "*                        step03                                      *"
# echo "**********************************************************************"
# apptainer exec instance://pipeline-inst Rscript --vanilla script/03_extract_experimental_uncertainties.R $3
# 
# echo "**********************************************************************"
# echo "*                        step04                                      *"
# echo "**********************************************************************"
# apptainer exec instance://pipeline-inst Rscript --vanilla script/04_tune_experimental_uncertainties_new.R $3
# echo "**********************************************************************"
# echo "*                        step05                                      *"
# echo "**********************************************************************"
# apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script/05_create_reference_jacobian_mpi.R $3

#echo "**********************************************************************"
#echo "*                        step06                                      *"
#echo "**********************************************************************"
#apptainer exec instance://pipeline-inst Rscript --vanilla script/06_tune_endep_hyperpars_with_prior.R $3
#echo "**********************************************************************"
#echo "*                        step07                                      *"
#echo "**********************************************************************"
##apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script/07_tune_talyspars.R $3
#apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script/07_tune_talyspars_all_endep_free.R $3
#echo "**********************************************************************"
#echo "*                        step07.5                                    *"
#echo "**********************************************************************"
#apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script/07_5_addGPobs.R $3
#echo "**********************************************************************"
#echo "*                        step08                                      *"
#echo "**********************************************************************"
###apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script/08_calculate_posterior_approximation_with_defect.R $3
#apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script/08_calculate_posterior_approximation.R $3
#echo "**********************************************************************"
#echo "*                        step09                                      *"
#echo "**********************************************************************"
#apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script/09_create_randomfiles_mpi.R $3
#echo "**********************************************************************"
#echo "*                        step10                                      *"
#echo "**********************************************************************"
##apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script/10_tune_talyspars_with_defect.R $3
#apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script/10_tune_talyspars_with_defect_gp_obs_all_endep_free.R $3
#echo "**********************************************************************"
#echo "*                        step11                                      *"
#echo "**********************************************************************"
#apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script/11_calculate_posterior_approximation_with_defect.R $3

echo "**********************************************************************"
echo "*                        step12                                      *"
echo "**********************************************************************"
#apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script/12_create_randomfiles.R $3
#apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script/12_create_randomfiles_new.R $3
apptainer exec instance://pipeline-inst mpirun -np 1 Rscript --vanilla script/12_create_randomfiles_test_energy-grid.R $3

apptainer instance stop pipeline-inst