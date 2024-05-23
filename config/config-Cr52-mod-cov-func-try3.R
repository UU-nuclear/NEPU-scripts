##################################################
#
#       CONFIGURATION OF PIPELINE
#
##################################################

# add a user library where we can install additional packages
# userLib <- "/TMC/alf/pipeline/eval-Cr-isotopes/eval-fe56-scripts/R-libs-user"
# .libPaths( c( .libPaths(), userLib) )

# working directory
workdir <- "/proj/naiss2024-22-324/ND-eval-pipeline/eval-fe56-scripts"
setwd(workdir)

source("config/required_packages.R")
source("config/required_sourcefiles.R")

tmp_dir <- file.path("/dev/shm",Sys.getenv("SLURM_JOB_ID"))

createTalysHandlers <- function() {

    # Initialize the talysR mpi interface

    # Important note: 1) Scripts that will run talys under mpi should call this function at the very
    #                    begining of the script. The execution of this function, blocks any futher
    #                    execution of the R-scripts on the slave ranks (rank-id>0). The R code in the
    #                    scripts should be executed only on the main thread, so call this as early
    #                    as possible!
    #                 2) TMPDIR = "/dev/shm" is an important specification because /dev/shm usually
    #                    resides in main memory. TALYS produces many thousand files per run
    #                    and normal disks and shared file systems cannot deal with this load
    #                    so it is a good idea to store them in main memory.
    #                 3) maxNumCPU set the number of requested talys workers. The number is an upper
    #                    limit on the number of workers. If maxNumCPU=0 the number of workers will be
    #                    the number of availible workers as given by the MPI interface.
    runOpts <- list(TMPDIR = tmp_dir)
    talysHnd <- initTALYSmpi(runOpts = runOpts, maxNumCPU=0)

    # initialize an alternative TALYS handler
    talysOptHnd <- createTalysFun(talysHnd, TMPDIR=tmp_dir)

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

# specify the reaction(s) to extract data fromthe EXFOR data base
# target reaction strings matching this regular expression
reacPat <- "\\(24-CR-52\\(N,[^)]+\\)[^,]*,,SIG\\)"

# should pipeline be executed with a very 
# small number of adjustable model parameters
# for testing purposes
fewParameterTest <- FALSE

# time interval to check for completed
# TALYS calculation in seconds
pollTime <- 1

# settings to retrieve EXFOR entries from
# the MongoDb database
mongo_dbname <- "exfor"
mongo_colname <- "entries"

# only use experimental data in that energy range
minExpEn <- 1.0
maxExpEn <- 50

# Specify energy grid for the final random files created in step 9.
# The grid used during the fit is based on this one, but limited to the range of
# the experimental data. Therefore the following must hold
# energyGridrandomFiles[1] < minExpEn
# energyGridrandomFiles[length(energyGridrandomFiles)] > maxExpEn
# Note that for proper error propagation the energy dependent paramters should cover the same range.
energyGridrandomFiles <- c(seq(0.1,1.,0.1),seq(1.2,8.0,0.2),seq(8.5,15.0,0.5),seq(16.,30,1),seq(32,60,2),seq(65,80,5),seq(90,200,10))

# moved the creation of the energy grid used in the fit to script 02, in order to limit it to where there is data


# default threshold energy for reaction channels
# if automatic determination fails
# (because of vanishing reaction cross section at all energies)
defaultThresEn <- 1 

# energy grid for energy-dependent TALYS parameters
# copy energyGridrandomFiles 
energyGridForParams <- energyGridrandomFiles
# take every 3rd point for the parameters
energyGridForParams <- energyGridForParams[seq(3,length(energyGridForParams),by=3)]
# limit at the maximum energy, including the energy grid-point above maxExpEn
energyGridForParams <- energyGridForParams[1:(max(which(energyGridForParams < maxExpEn))+1)]
# add the energy at "zero" (must be a positive energy, otherwise talys ignores it)
energyGridForParams <- c(1.e-06,energyGridForParams)

# specify which paramters to make energy dependent
# a data.frame object named enParDt containing the columns par and proj
# is used for this, for example
# > enParDt
#           par proj
# 1    v1adjust    n
# 2    d1adjust    n
tmpPar <- paste0(c('v1','d1','w1','vso1','wso1','rc','av','avd','avso','aw','awso','rv','rvd','rvso','rwd','rwso','awd'),'adjust')
tmpProj <- c('n','p','d','t','h','a')
enParDt <- data.table(expand.grid(par = tmpPar, proj = tmpProj))
#enParDt <- enParDt[!(par=='rcadjust' & proj=='n')] # remove Coloumb radius for the neutron

# specification of the TALYS input file used as template
# param_template_path <- file.path(workdir,"indata/n_Fe_056.inp")
# I think that it should not really matter that this file is for Fe-56, only the parameters are extracted from the file
# the target and projectile are specified sepparately
# the input will be searched for in the indata directory, if not found there, it will be downloaded from
# https://tendl.web.psi.ch/tendl_2019/
# use the following keywords to specify which nuclide and projectile
tendl_element <- "Cr"
tendl_mass <- 52
tendl_projectile <- "n"


# optional argument to set the finite difference used by talys to calculate the Jacobian
# default value is talys_finite_diff <- 0.01
talys_finite_diff <- 0.01

# random generator seed for optimization of experimental uncertainties
# impacts the initial extra uncertainties in the optimization setup
tuneExpUncSeed <- 11

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
#outdataPath <- file.path(workdir, "/outdata-try1")
outdataPath <- file.path(workdir, "/outdata-mod-cov-func-try3")
dir.create(outdataPath, recursive=TRUE, showWarnings=FALSE)

# specify the directory were status information and plots during the 
# optimization using the Levenberg-Marquardt algorithm should be stored
savePathLM <- file.path(outdataPath, "/LMalgo")

# random seed to create TALYS randomfiles
talysFilesSeed <- 13

# number of TALYS randomfiles to be created
 numTalysFiles <- 1000

# where to store the TALYS results
# content of TALYS result directories is stored as tar archives
# needed for the creation of ENdf randomfiles using modified TASMAN
pathTalys <- file.path(outdataPath, "random-files")
savePathTalys <- pathTalys

# where to save plots produced by the scripts in eval-fe56/script/visualization
plotPath <- file.path(outdataPath, '/plots')
