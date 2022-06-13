##################################################
#
#       CONFIGURATION OF PIPELINE
#
##################################################

# working directory
workdir <- "/home/alfgo462/NucDat/pipeline/eval-fe56-singularity/workdir/"
setwd(workdir)

source("/opt/pipeline/eval-fe56/required_packages.R")
source("/opt/pipeline/eval-fe56/required_sourcefiles.R")
library(clusterTALYSmpi);

# temporarily place the sourcing of the clusterTALYSmpi script here
# later on this should be installed as a package
#source("/home/alf/programs/clusterTALYSmpi/R/clusterTALYSmpi.R")

# should pipeline be executed with a very 
# small number of adjustable model parameters
# for testing purposes
fewParameterTest <- TRUE
max_number_of_free_pars <- 20

# should TALYS calculations be performed
# within the Docker container
containerTest <- FALSE

# a valid path on the remote machine and the
# local machine, respectively, to store
# result files 
calcdir_rem <- "/TMC/alf/eval-fe56-singularity/remCalcDir"
calcdir_loc <- paste0(workdir,"/localCalcDir")

# time interval to check for completed
# TALYS calculation in seconds
pollTime <- 10

# settings to retrieve EXFOR entries from
# the MongoDb database
mongo_dbname <- "exfor"
mongo_colname <- "entries"

# energy grid for TALYS calculations
energyGrid <- seq(0.1, 30.001, length = 50)
#energyGrid <- c(seq(0.1,10,by=0.3),seq(10.6,20,by=0.6),seq(20.2,30.001,length=11))

# Option to specify different grid for the final random files created i step 9. 
# For example an extrapolation from data region.
# Note that for proper error propagation the energy dependent paramters should cover the same range.
# energyGridrandomFiles <- seq(0.1, 30.001, length = 100)
energyGridrandomFiles <- energyGrid

# default threshold energy for reaction channels
# if automatic determination fails
# (because of vanishing reaction cross section at all energies)
defaultThresEn <- 1 

# energy grid for energy-dependent TALYS parameters
energyGridForParams <- seq(0,31,by=2)

# specification of the TALYS input file used as template
# param_template_path <- "/opt/pipeline/eval-fe56/indata/n_Fe_056.inp"
param_template_path <- "/home/alfgo462/NucDat/pipeline/eval-fe56-singularity/workdir/indata/n_Fe_056.inp"

# instantiate the transformation used for all parameters of the form ...adjust
# parameters are restricted to the interval (0.5, 1.5), in other words:
# the maximal deviation from the default values is 50%
paramTrafo <- generateTrafo(1, 0.5, 4) 

# random generator seed for optimization of experimental uncertainties
# impacts the initial extra uncertainties in the optimization setup
tuneExpUncSeed <- 11

# only use experimental data in that energy range
minExpEn <- 2
maxExpEn <- 30

# set up the handlers to map TALYS results to EXFOR entries
subentHandler <- createSubentHandler(createDefaultSubentHandlerList())
exforHandler <- createExforHandler(subentHandler)
# abuAgent <- createAbuAgent("talys/structure/abundance/")
# subentHandler$getHandlerByName("handler_ntot_nat")$configure(list(abuAgent = abuAgent))

# maximum number of iterations for Levenberg-Marquardt algorithm
maxitLM <- 30

# if the relative difference between subsequent iterations of the
# Levenberg-Marquardt algorithm falls below this value,
# the optimization procedure terminates
reltolLM <- 1e-5

# where to save output data
outdataPath <- file.path(workdir, "/outdata_testLM")
#outdataPath <- file.path(workdir, "/outdata_new_grid")
dir.create(outdataPath, recursive=TRUE, showWarnings=FALSE)

# specify the directory were status information and plots during the 
# optimization using the Levenberg-Marquardt algorithm should be stored
savePathLM <- file.path(outdataPath, "/LMalgo")

# random seed to create TALYS randomfiles
talysFilesSeed <- 13

# number of TALYS randomfiles to be created
 numTalysFiles <- 300
 #numTalysFiles <- 10

# where to store the TALYS results on the remote machine
# content of TALYS result directories is stored as tar archives
# needed for the creation of ENdf randomfiles using modified TASMAN
# pathTalys <-paste0(workdir,"/talysResults")
pathTalys <-paste0("/tmp/talysResults")
# pathTalys <- "/TMC/alf/eval-fe56-singularity/talysResults/full-eval"
savePathTalys <- pathTalys

# where to save plots produced by the scripts in eval-fe56/script/visualization
plotPath <- file.path(outdataPath, '/plots')


createTalysHandlers <- function() {

    # Important note: 1) TMPDIR = "/dev/shm" is an important specification because /dev/shm usually
    #                    resides in main memory. TALYS produces many thousand files per run
    #                    and normal disks and shared file systems cannot deal with this load
    #                    so it is a good idea to store them in main memory.
    #                 2) bindir points to the location of the runTALYSmpi program that takes care 
    #                    of executing talys wiht the mpi-interface.
    #                     i)   if the mpi session is invoked inside the container, for example
    #                             singularity exec <sif-file> mpirun -np 1 Rscript --vanilla <script> <config>
    #                          it is not necesary to provide this, since the code will find it in
    #                          /usr/local/bin
    #                     ii)  if the mpi session is invoked outside the container, for example
    #                             mpirun -np 1 exec <sif-file> Rscript --vanilla <script> <config>
    #                          runTALYSmpi must be compiled outside the container, and bindir should be 
    #                          pointing to the location of the reulting binary
    #                 3) maxNumCPU set the number of requested talys workers. The number is an upper
    #                    limit on the number of workers. If maxNumCPU=0 the number of workers will be
    #                    the number of availible workers as given by the MPI interface.
    runOpts <- list(TMPDIR = "/dev/shm/talysTemp",bindir = "/usr/local/bin")
    talysHnd <- initClusterTALYSmpi(talysExe = "talys", runOpts = runOpts, maxNumCPU=0)

    # initialize an alternative TALYS handler
    talysOptHnd <- createTalysFun(talysHnd)

    # Difference between talysHnd and talysOptHnd:
    #   talysHnd is a lower-level interface that provides
    #            the functions run, isRunning, and result.
    #            The input specification is passed as a list
    #            with input keywords and values and the output
    #            specification as a datatable enumerating the
    #            observables of interest

    #   talysOptHnd provides the functions fun and jac which 
    #               take a vector x as input and return either
    #               a vector of observables (fun) or the Jacobian 
    #               matrix (jac). Default parameter values and
    #               which values are present in x is specified
    #               via additional setter functions. Functions
    #               provided by talysOptHnd rely on those 
    #               provided by talysHnd.

    list(talysHnd = talysHnd,
         talysOptHnd = talysOptHnd)
}