library(data.table)
library(talysR)
#library(doMPI)

runOpts <- list()
runOpts$TMPDIR <- "/dev/shm/talysTemp"
talysClustObj <- initTALYSmpi(maxNumCPU=0, runOpts=runOpts)

energyGrid <- c(seq(0.1,1.,0.1),seq(1.2,8.0,0.2),seq(8.5,15.0,0.5),seq(16.,30,1),seq(32,40,2))

paramList <- list(projectile = "n",
                   element = "Cr",
                   mass = 52,
                   energy = energyGrid)

outSpec <- data.table(REAC = "CS/TOT",
                      L1 = seq(1, 3, length = 10),
                      L2 = 0, L3 = 0) 

startTime <- Sys.time()
runObj <- talysClustObj$run(paramList, outSpec, saveDir="/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/defect-model")
stopTime <- Sys.time()

exec_time <- as.double(stopTime-startTime,units="mins")

cat("talys execution time: ",exec_time,"\n")

cat("the end\n")
