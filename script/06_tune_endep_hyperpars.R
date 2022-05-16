#
# DESCRIPTION OF STEP
#
# tune the hyperparameters of the Gaussian
# processes attached to the energy-dependent
# TALYS parameters using maximunm likelihood
# optimization.
#

#################################################
#       SCRIPT Setup
##################################################

args = commandArgs(trailingOnly=TRUE)


if (length(args)==0) {
  source("./config/config.R")
  stop("No config file supplied, using default file config.R", call.=FALSE)
} else if (length(args) > 1) {
  stop("Script only accepts one argument.", call.=FALSE)
} else {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}


#################################################
#       SCRIPT PARAMETERS
##################################################

scriptnr <- 6L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

subents <- read_object(1, "subents")
refParamDt <- read_object(2, "refParamDt")
extNeedsDt <- read_object(2, "extNeedsDt")
modList <- read_object(3, "modList")
expDt <- read_object(3, "expDt")
updSysDt <- read_object(4, "updSysDt")
fullSensDt <- read_object(5, "fullSensDt") 

##################################################
#       START OF SCRIPT
##################################################

# define objects to be returned
outputObjectNames <- c("optExpDt", "optSysDt", "optGpDt")
check_output_objects(scriptnr, outputObjectNames)

# prepare the data table with systematic components
# this comprises both systematic experimental 
# errors and the model parameters

talysHandler <- createSysCompModelHandler()
talysHandler$setRef(extNeedsDt, fullSensDt, refParamDt,
                    exforHandler, c(subents, modList$SUBENT))
talysHandler$setPrior(refParamDt)
modelSysDt <- talysHandler$createSysDt()

setkey(updSysDt, IDX)
curSysDt <- rbind(updSysDt, modelSysDt, fill=TRUE)
curSysDt <- curSysDt[!grepl("REACEXP-", EXPID)]
curSysDt[, IDX := seq_len(.N)]
curSysDt[, ADJUSTABLE := FALSE]

# set up Gaussian processes for energy dependent parameters
gpHandler <- createSysCompGPHandler()
parnames_endep <- unique(curSysDt[ERRTYPE == "talyspar_endep", EXPID])
for (curParname in parnames_endep) {
    gpHandler$addGP(curParname, 0.1, 2, 1e-4)
}

# create a mapping matrix from the model output
# to the grid of the experimental data
Sexp <- exforHandler$getJac(expDt, extNeedsDt, subents) 

# prepare the experimental data
curExpDt <- copy(expDt)
setkey(curExpDt, IDX)
setkey(extNeedsDt, IDX)
# REFDATA contains the cross sections of the reference model
# calculation. It serves as the expansion point in the
# Taylor approximation of the model entering the Bayesian
# update
curExpDt[, REFDATA := as.vector(Sexp %*% extNeedsDt$V1)]
curExpDt[, ADJUSTABLE := FALSE]

gpDt <- gpHandler$createGPDt()
gpHandler$updateSysDt(curSysDt)

# define which GP hyperparameters are adjustable
gpDt[, ADJUSTABLE := PARNAME %in% c("sigma","len")]
gpDt[, IDX := seq_len(.N)]

# initialize handler for systematic experimental uncertainties
# TODO: It should not be needed to define the systematic error
#       components during the setup. The category (column name)
#       should be enough.
normHandler <- createSysCompNormHandler("DATAREF")
tmpDt <- updSysDt[grepl("EXPID-", EXPID)]
expIds <- sub("EXPID-", "", tmpDt$EXPID)
normHandler$addSysUnc("EXPID", expIds, 0, tmpDt$UNC,
                      rel = (tmpDt$ERRTYPE == "sys-rel"))

# create global handler and register the individual handlers
sysCompHandler <- createSysCompHandler()
sysCompHandler$addHandler(normHandler)
sysCompHandler$addHandler(talysHandler)
sysCompHandler$addGPHandler(gpHandler)

# setup optimization specification

optfuns <- createMLOptimFuns()
optfuns$setDts(curExpDt, curSysDt, gpDt,
               sysCompHandler = sysCompHandler)

# define the hyperparameter limits
setkey(gpDt, IDX)
lowerLims <- rep(NA_real_, nrow(gpDt))
upperLims <- rep(NA_real_, nrow(gpDt))
lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- 0.1
lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- 2
lowerLims <- lowerLims[gpDt$ADJUSTABLE]
upperLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- 0.5
upperLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- 50
upperLims <- upperLims[gpDt$ADJUSTABLE]

# initial configuration
#initPars <- lowerLims + runif(length(lowerLims)) * abs(upperLims - lowerLims)
initPars <- lowerLims

# optimize hyperparameters
# this may take a few minutes
optRes <- optim(initPars, optfuns$logLike, optfuns$gradLogLike, method = "L-BFGS-B",
                lower = lowerLims, upper = upperLims, control = list(fnscale = -1))

newDts <- optfuns$getModifiedDts(optRes$par)

optExpDt <- newDts$expDt
optSysDt <- newDts$sysDt
optGpDt <- newDts$gpDt

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)

