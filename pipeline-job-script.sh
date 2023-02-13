#!/bin/bash -l

# this job script is specifical for step 02 of the pipeline that does not use
# more than two cores:
# one for the main thread and one worker to do the talys calculation
# so I will run it on a single node requsting two cores

#SBATCH -A naiss2023-22-584
#SBATCH -p core -n 2
#SBATCH -t 40:00
#SBATCH -J talysTemp

#don't I need to load the singularity module?
apptainer exec /proj/naiss2023-22-58/ND-eval-pipeline mpirun -np 1 Rscript --vanilla script-Cr/02_create_reference_calculation.R /proj/naiss2023-22-58/ND-eval-pipeline/eval-fe56-scripts/config/config-test.R