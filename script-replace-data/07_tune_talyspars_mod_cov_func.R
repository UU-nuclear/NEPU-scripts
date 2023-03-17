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

talysHnds <- createTalysHandlers()

library(stringr)


#################################################
#       SCRIPT PARAMETERS
##################################################

scriptnr <- 7L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

subents <- read_object(4, "fake_subents")
refParamDt <- read_object(2, "refParamDt")
#extNeedsDt <- read_object(2, "extNeedsDt")
extNeedsDt <- read_object(4, "extNeedsDt")
modList <- read_object(3, "modList")
full_covMat <- read_object(4, "full_covMat")
fullSensDt <- read_object(5, "fullSensDt") 
optParamDt <- read_object(5, "optParamDt")
mask <-  read_object(5,"mask")
Sexp <-  read_object(5,"Sexp")
optExpDt <- read_object(6, "optExpDt")
optSysDt <- read_object(6, "optSysDt")
optGpDt <- read_object(6, "optGpDt")

extNeedsDt[,IDX:=seq_len(.N)]

##################################################
#       START OF SCRIPT
##################################################
print("-----------------------------------------------------")
print("----------------------script 07----------------------")
print("-----------------------------------------------------")

# define objects to be returned
outputObjectNames <- c("optRes", "refPar", "P0", "yexp", "D", "S0", "X",
                       "optSysDt_allpars", "optSysDt_optpars")
check_output_objects(scriptnr, outputObjectNames)

talys <- talysHnds$talysOptHnd

# create the parameter transformation object
# set the parameter transformation to be centered at the prior mean/mode
# the ranges of the parameters are the ones specified in the TALYS manual
paramTrafo <- parameterTransform(
                  x0 = unlist(optParamDt[ADJUSTABLE==TRUE,PARVAL]),
                  delta = optParamDt[ADJUSTABLE==TRUE,unlist(PARVAL) - PARMIN])

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

# sanity check: check if we have set up the model correctly
# refPar <- optParamDt[ADJUSTABLE == TRUE, unlist(PARVAL)]
# fn <- talys$fun(refPar)

# we remove ADJUST columns not needed in this step.
# The ADJUST column in optParamDt defines
# which TALYS parameters (both energy-dependent and 
# -independent) should be optimized
optExpDt[, ADJUSTABLE := NULL]
optSysDt[, ADJUSTABLE := NULL]
optGpDt[, ADJUSTABLE := NULL]

# Both optParamDt and optSysDt contain the TALYS parameters.
# We need the parameters in the format of optSysDt to 
# calculate the prior covariance matrix for the model 
# parameters. However, for the optimization using the LM
# algorithm we need them in the form of optParamDt.
# As an additional difficulty, the parameter specification
# in the column EXPID in optSysDt is different from that
# in column PARNAME in optParamDt for energy-dependent parameters.
# Here is some ad-hoc code to map between.
# TODO: Improvements for the future
#       sysCompTalysHandler should be changed to be compatible
#       with parameter specifications as in the form of optParamDt
#       Conversion between 'numeric' and 'character' may
#       cause problems in identifying corresponding components.
optSysDt[grepl("^TALYS-", EXPID), PARNAME := {
    parname <- sub("^TALYS-","", EXPID)
    ifelse(is.na(EN), parname, sub("adjust", paste0("adjust(", energyGridForParams[EN], ")"), parname))
}, by="IDX"]
setkey(optSysDt, PARNAME)
setkey(optParamDt, PARNAME)
tmpIdcs <- match(optParamDt[ADJUSTABLE==TRUE,PARNAME], optSysDt$PARNAME)
stopifnot(!is.unsorted(tmpIdcs))
stopifnot(length(tmpIdcs) == nrow(optParamDt[ADJUSTABLE==TRUE]))

# keep a reduced version of optSysDt
# with TALYS parameters that were initially considered
# for optimization
optSysDt_allpars <- optSysDt[grepl("^TALYS-", EXPID),]
setkey(optSysDt_allpars, IDX)
optSysDt_allpars[, IDX:=seq_len(.N)]

# remove TALYS parameters from optSysDt that will not be optimized
# note that optSysDt still also contains systematic experimental components
fixedParnames <- optParamDt[ADJUSTABLE == FALSE, PARNAME]
optSysDt <- optSysDt[! PARNAME %in% fixedParnames]

# store the TALYS parameters to be optimized in
# an additional data table
optSysDt_optpars <- optSysDt_allpars[! PARNAME %in% fixedParnames]

# recreate the index because we have deleted rows
setkey(optSysDt, IDX)
setkey(optParamDt, IDX)
optSysDt[, IDX := seq_len(.N)]
stopifnot(optSysDt[grepl("TALYS",EXPID),PARNAME] == optParamDt$PARNAME[optParamDt$ADJUSTABLE == TRUE])

# prepare the Gaussian process handler 
gpHandler <- createSysCompGPHandler()

# prepare the handlers to map systematic uncertainties of the experiments
#normHandler <- createSysCompNormHandler("DATAREF")
#normHandler$addSysUnc("EXPID", "", 0, 0, FALSE)

# prepare the TALYS handler to map from model parameters to predictions
talysHandler <- createSysCompModelHandler()
talysHandler$setRef(extNeedsDt, fullSensDt, refParamDt,
                    exforHandler, c(subents, modList$SUBENT))
talysHandler$setPrior(refParamDt)


# create global handler and register the individual handlers
sysCompHandler <- createSysCompHandler()
#sysCompHandler$addHandler(normHandler)
sysCompHandler$addHandler(talysHandler)
sysCompHandler$addGPHandler(gpHandler)

# construct the matrix to map from systematic error components
# of the experiments to the measurement points
S <- sysCompHandler$map(optExpDt, optSysDt, ret.mat = TRUE)
# construct the covariance matrix containing both the 
# covariances of the systematic experimental errors and
# the covariances for the model parameters. This matrix is diagonal
# for the portions relating to experiments and energy-independent
# parameters. The blocks associated with energy-dependent TALYS
# parameters are constructed by a Gaussian process and contain
# therefore correlations.
P <- sysCompHandler$cov(optSysDt, optGpDt, ret.mat = TRUE)

# The Levenberg-Marquardt routine assumes that systematic
# components in optSysDt are only related to experiments
# and not to model parameters. The latter were introduced
# in step 06 to optimize the hyperparameters of the GPs
# associated with energy-dependent TALYS parameters.
# Therefore we need to remove them now.
setkey(optSysDt, IDX)
setkey(optExpDt, IDX)
expSel <- optSysDt[, !grepl("TALYS-", EXPID)]
talysSel <- !expSel

P0 <- P[talysSel, talysSel]
X <- P[expSel, expSel] 
S0 <- S[, expSel]
# hack to not map the systematic experimental errors, which are already taken into account
# in full_covMat: set all elemtents to 0
S0[,] <- 0

#D <- Diagonal(x = optExpDt$UNC^2)
D <- full_covMat
yexp <- getDt_DATA(optExpDt)
setkey(optParamDt, IDX)
refPar <- optParamDt[ADJUSTABLE==TRUE, unlist(PARVAL)]

# because we removed rows from optSysDt, we need to restore
# a continuous index which is required for functions that create
# matrices and rely on a continuous and complete set of
# indices
#optSysDt <- optSysDt[sel,]
#setkey(optSysDt, IDX)
#optSysDt[, IDX := seq_len(.N)]

if (!dir.exists(savePathLM)) dir.create(savePathLM, recursive=TRUE)
#loggerLM <- createLoggerLM(talys, savePathLM)
loggerLM <- createLoggerLMalt(savePathLM)

# uncomment the line below to start from last parameterset of previous LM run
#pinit <- read_object(7, "optRes")$par
pinit <- refPar

# The full Jacobian was calculated in step 05, no need to recalculate it
Spar <- with(fullSensDt,
             sparseMatrix(i = IDX1, j = IDX2, x = X,
                          dims = c(nrow(extNeedsDt), nrow(refParamDt))))
#Sexp <- exforHandler$getJac(expDt, extNeedsDt, subents)
Sglob <- Sexp %*% Spar 
Jinit <- Sglob[,optParamDt[ADJUSTABLE==TRUE]$IDX] # should now be read as Sopt from select_parameters.R

#cat("Started calculations at", as.character(Sys.time()), "\n")  
#optRes <- LMalgo(talys$fun, talys$jac, pinit = pinit, p0 = refPar, P0 = P0, D = D, S = S0, X = X, yexp =yexp,
#                 lower = rep(-Inf, length(refPar)), upper = rep(Inf, length(refPar)), logger = loggerLM,
#                 control = list(maxit = maxitLM, reltol = reltolLM, acc = FALSE, alpha=0.75, acc_step = 1e-1))
#cat("Finished calculations at", as.character(Sys.time()), "\n")

startTime <- Sys.time()
cat("Started calculations at", as.character(startTime), "\n")
source("LMalgo_parallel/LMalgo_parallel.R")
optRes <- LMalgo_parallel(talys$fun, talys$jac, pinit = pinit, p0 = refPar, P0 = P0, D = D, S = S0, X = X, yexp =yexp,
                 lower = rep(-Inf, length(refPar)), upper = rep(Inf, length(refPar)), logger = loggerLM,
                 control = list(maxit = maxitLM, reltol = reltolLM, steptol=0.1*talys$getEps(), acc = FALSE, alpha=0.75, acc_step = 1e-1, nproc = 31, strategy = "gain"),J=Jinit)

stopTime <- Sys.time()
cat("Finished calculations at", as.character(stopTime), "\n")

exec_time <- as.double(stopTime-startTime,units="hours")

cat("total execution time: ",exec_time," hours\n")

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)
