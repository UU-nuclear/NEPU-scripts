#
# DESCRIPTION OF STEP
#
# Use the covariance matrices constructed
# during the previous steps and solve the
# non-linear Generalized-Least-Squares problem
# with the Levenberg-Marquardt algorithm.
#

#################################################
#       SCRIPT Setup
##################################################

args = commandArgs(trailingOnly=TRUE)



if(length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
} else if (length(args) > 1) {
  stop("Script only accepts one argument.", call.=FALSE)
}



#################################################
#       SCRIPT PARAMETERS
##################################################

scriptnr <- 10L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

subents <- read_object(1, "subents")
refParamDt <- read_object(2, "refParamDt")
extNeedsDt <- read_object(2, "extNeedsDt")
modList <- read_object(3, "modList")
fullSensDt <- read_object(5, "fullSensDt") 
optExpDt <- read_object(6, "optExpDt")
optSysDt <- read_object(6, "optSysDt")
optGpDt <- read_object(6, "optGpDt")
P0 <- read_object(7, "P0")
yexp <- read_object(7, "yexp")
D <- read_object(7, "D")
S0 <- read_object(7, "S0k")
X <- read_object(7, "Xk")
optParamDt <- read_object(7, "optParamDt")
mask <- read_object(7, "mask")
Sexp <- read_object(7, "Sexp")
optSysDt_allpars <- read_object(7,"optSysDt_allpars")
optSysDt_optpars <- read_object(7,"optSysDt_optpars")
#finalPars <- read_object(8, "finalPars")
refPar <- read_object(7,"refPar")
##################################################
#       START OF SCRIPT
##################################################

# TODO: Right now this script is basically a copy paste of step07.
# Should implement that the results from step07 about which parameters to keep adjustable is read
# from step07 instead of redoing the calculation. This is mainly for cleanliness but also makes 
# the execution less prone to errors since if I want to change something I don't need to do it in two 
# sepparate scripts.

# define objects to be returned
outputObjectNames <- c("optRes", "optParamDt", "Sexp", "mask",
                       "refPar", "P0", "yexp", "D", "S0", "X",
                       "optSysDt_allpars", "optSysDt_optpars")
check_output_objects(scriptnr, outputObjectNames)

# set the parameter transformation to be centered at the prior mean/mode
# the ranges of the parameters are the ones specified in the TALYS manual
paramTrafo <- parameterTransform(
                  x0 = unlist(optParamDt[ADJUSTABLE==TRUE,PARVAL]),
                  delta = optParamDt[ADJUSTABLE==TRUE,unlist(PARVAL) - PARMIN])

talysHnds <- createTalysHandlers()
talys <- talysHnds$talysOptHnd

# define the default parameters
# the ADJUST column in optParamDt determines whether the 
# associated parameter is considered adjustable or fixed
talys$setPars(optParamDt)
# the optimization will be undertaken using transformed
# TALYS parameters to make sure that parameters remain
# within given limits. paramTrafo is defined in config.R
talys$setParTrafo(paramTrafo$fun, paramTrafo$jac)
# define the required predictions from TALYS in order to
# map to the grid of the experimental data given in EXFOR
talys$setNeeds(extNeedsDt)
# define the matrix to map from the output observables of
# TALYS to the information given in EXFOR
talys$setSexp(Sexp)
# knowing that some parameters have not any impact on 
# the comparison between model prediction and experiment
# we exclude the calculation of the associated elements
# in the Jacobian matrix.
talys$setMask(mask)
# the finite difference used to calculate numerically the
# Jacobian matrix
if(!exists("talys_finite_diff")) talys_finite_diff <- 0.01
talys$setEps(talys_finite_diff)

if (!dir.exists(savePathLM)) dir.create(savePathLM, recursive=TRUE)
loggerLM <- createLoggerLMalt(savePathLM)

# uncomment the line below to start from last parameterset of previous LM run
optRes7 <- read_object(7, "optRes")
pinit <- optRes7$par
Jinit <- optRes7$jac

cat("Started calculations at", as.character(Sys.time()), "\n")  
source("LMalgo_parallel/LMalgo_parallel.R")
optRes <- LMalgo_parallel(talys$fun, talys$jac, pinit = pinit, p0 = refPar, P0 = P0, D = D, S = S0, X = X, yexp =yexp,
                 lower = rep(-Inf, length(refPar)), upper = rep(Inf, length(refPar)), logger = loggerLM,
                 control = list(maxit = maxitLM, reltol = reltolLM, steptol=0.1*talys$getEps(), acc = FALSE,
                  alpha=0.75, acc_step = 1e-1, nproc = 29, strategy = "gain", mu=10),J=Jinit)
cat("Finished calculations at", as.character(Sys.time()), "\n")

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)

