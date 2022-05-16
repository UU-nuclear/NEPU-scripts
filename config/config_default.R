##################################################
#
#       CONFIGURATION OF PIPELINE
#
##################################################

# root directory of the evaluation pipeline
rootpath <- "/home/user/eval-fe56"
setwd(rootpath)

source("required_packages.R")
source("required_sourcefiles.R")

# should pipeline be executed with a very 
# small number of adjustable model parameters
fewParameterTest <- FALSE

# should TALYS calculations be performed
# within the Docker container
containerTest <- FALSE

# number of CPUs limiting the number of 
# TALYS calculations performed in parallel
maxNumCPU <- as.numeric(Sys.getenv('maxNumCPU'))
if (is.na(maxNumCPU)) maxNumCPU <- 32

# settings needed to connect to the cluster
# and run TALYS calculations in parallel
ssh_login <- "user@host"
ssh_pw <- "password"
cluster_node <- "host"

# a valid path on the remote machine and the
# local machine, respectively, to store
# result files 
calcdir_rem <- "/home/user/remcalcdir"
calcdir_loc <- "/home/user/calcdir"

# time interval to check for completed
# TALYS calculation in seconds
pollTime <- 60

# settings to retrieve EXFOR entries from
# the MongoDb database
mongo_dbname <- "exfor"
mongo_colname <- "entries"

# energy grid for TALYS calculations
energyGrid <- seq(0.1, 30.001, length = 100)

# Option to specify different grid for the final random files created i step 9. 
# For example an extrapolation from data region.
# Note that for proper error propagation the energy dependent paramters should cover the same range.
energyGridrandomFiles <- energyGrid

# default threshold energy for reaction channels
# if automatic determination fails
# (because of vanishing reaction cross section at all energies)
defaultThresEn <- NA

# energy grid for energy-dependent TALYS parameters
energyGridForParams <- seq(0,31,by=2)

# specification of the TALYS input file used as template
param_template_path <- file.path(rootpath, "indata/n_Fe_056.inp")

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
outdataPath <- file.path(rootpath, "outdata_test")
dir.create(outdataPath, recursive=TRUE, showWarnings=FALSE)

# specify the directory were status information and plots during the 
# optimization using the Levenberg-Marquardt algorithm should be stored
savePathLM <- file.path(rootpath, "outdata/log/LMalgo")

# random seed to create TALYS randomfiles
talysFilesSeed <- 13

# number of TALYS randomfiles to be created
numTalysFiles <- 300

# where to store the TALYS results on the remote machine
# content of TALYS result directories is stored as tar archives
# needed for the creation of ENdf randomfiles using modified TASMAN
pathTalys <-"/home/user/talysResults"

savePathTalys <-paste0(cluster_node,":",pathTalys)
# where to save plots produced by the scripts in eval-fe56/script/visualization
plotPath <- file.path(outdataPath, 'plots')


createTalysHandlers <- function() {
  
  # wait to avoid some problem due to connecting
  # too quickly in sequence
  Sys.sleep(1)
  
  # set up the connection to the cluster
  # and the functionality to run TALYS in parallel
  remHnd <- initSSH(ssh_login, ssh_pw,
                    tempdir.loc = calcdir_loc,
                    tempdir.rem = calcdir_rem)
  
  clustHnd <- initCluster(functions_multinode, remFun = remHnd) 
  
  # Important note: TMPDIR = "/dev/shm" is an important specification because /dev/shm usually
  #                 resides in main memory. TALYS produces many thousand files per run
  #                 and normal disks and shared file systems cannot deal with this load
  #                 so it is a good idea to store them in main memory.
  talysHnd <- initClusterTALYS(clustHnd, talysExe = "/usr/local/bin/talys", calcsPerJob = 1000,
                               runOpts = list(TMPDIR = "/dev/shm/talysTemp"))
  
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
  list(remHnd = remHnd,
       clustHnd = clustHnd,
       talysHnd = talysHnd,
       talysOptHnd = talysOptHnd)
}


# function that copies worker controller script to cluster
# and returns the command that needs to be run on worker nodes
# to watch for and execute transmitted jobs
slaveSetupCmd <- function(nohup = TRUE, launch = FALSE) {
  con <- createTalysHandlers()
  cmdstr <- con$clustHnd$startNodeController(con$clustHnd, maxNumCpus=maxNumCPU)
  con$clustHnd$closeCon()
  if (isTRUE(launch))
    system(cmdstr, wait = FALSE, ignore.stderr = TRUE)
  if (isTRUE(nohup))
    cmdstr <- paste0("nohup ", cmdstr, " &")
  cmdstr
}

if (isTRUE(containerTest)) {
  slaveSetupCmd(nohup = FALSE, launch = TRUE)
}

