#!/bin/bash -l

# submit a series of slurm jobs to perform the different steps in the pipeline, when one
# job depends on the results of another, make sure that they are executed in the correct
# sequence using the dependency flag of sbatch

JOBID_01=$(sbatch pipeline-job-01.sh) &&
JOBID_02=$(sbatch --dependency=afterok:${JOBID_01##* } pipeline-job-02.sh) &&
JOBID_03=$(sbatch --dependency=afterok:${JOBID_02##* } pipeline-job-03.sh) &&
JOBID_04=$(sbatch --dependency=afterok:${JOBID_03##* } pipeline-job-04.sh) &&
JOBID_05=$(sbatch --dependency=afterok:${JOBID_04##* } pipeline-job-05.sh) &&
JOBID_06=$(sbatch --dependency=afterok:${JOBID_05##* } pipeline-job-06.sh) &&
JOBID_07=$(sbatch --dependency=afterok:${JOBID_07##* } pipeline-job-07.sh)

echo $JOBID_01
echo $JOBID_02
echo $JOBID_03
echo $JOBID_04
echo $JOBID_05
echo $JOBID_06
echo $JOBID_07

#############################

# JOBID_04=$(sbatch pipeline-job-04.sh) &&
# JOBID_05=$(sbatch --dependency=afterok:${JOBID_04##* } pipeline-job-05.sh) &&
# JOBID_06=$(sbatch --dependency=afterok:${JOBID_05##* } pipeline-job-06.sh)
# 
# echo $JOBID_04
# echo $JOBID_05
# echo $JOBID_06

#############################
# 
# JOBID_03=$(sbatch pipeline-job-03.sh) &&
# JOBID_04=$(sbatch --dependency=afterok:${JOBID_03##* } pipeline-job-04.sh) &&
# JOBID_05=$(sbatch --dependency=afterok:${JOBID_04##* } pipeline-job-05.sh) &&
# JOBID_06=$(sbatch --dependency=afterok:${JOBID_05##* } pipeline-job-06.sh)
# 
# echo $JOBID_03
# echo $JOBID_04
# echo $JOBID_05
# echo $JOBID_06