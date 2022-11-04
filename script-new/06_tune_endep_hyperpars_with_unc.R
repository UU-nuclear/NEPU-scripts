#
# DESCRIPTION OF STEP
#
# tune the hyperparameters of the Gaussian
# processes attached to the energy-dependent
# TALYS parameters using maximunm likelihood
# optimization.
#
# This script is the same as 06_tune_endep_hyperpars_fixed_sigma.R
# with the difference that the upper limit for sigma is set by the prior
# uncertainty of the corresponding energy independent parameter

#################################################
#       SCRIPT Setup
##################################################

args = commandArgs(trailingOnly=TRUE)

if (length(args)==1) {
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
optParamDt <- read_object(6, "optParamDt")
extNeedsDt <- read_object(2, "extNeedsDt")
modList <- read_object(3, "modList")
expDt <- read_object(3, "expDt")
updSysDt <- read_object(4, "updSysDt")
fullSensDt <- read_object(6, "fullSensDt") 

##################################################
#       START OF SCRIPT
##################################################



print("-----------------------------------------------------")
print("----------------------script 06----------------------")
print("-----------------------------------------------------")

# define objects to be returned
outputObjectNames <- c("optExpDt", "optSysDt", "optGpDt","sysCompHandler")
check_output_objects(scriptnr, outputObjectNames)

# prepare the data table with systematic components
# this comprises both systematic experimental 
# errors and the model parameters

# For the talys handler to work we need to set all parameters in optParamDt
# that appear in FullSensDt to be adjustable
optParamDt[,ADJUSTABLE:=FALSE]
optParamDt[IDX %in% unique(fullSensDt[,IDX2]),ADJUSTABLE:=TRUE]

talysHandler <- createSysCompModelHandler()
talysHandler$setRef(extNeedsDt, fullSensDt, optParamDt,
                    exforHandler, c(subents, modList$SUBENT))
talysHandler$setPrior(optParamDt)
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
    gpHandler$addGP(curParname, 0.5*curSysDt[EXPID==curParname,UNC][1], 50, 1e-4)
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
gpDt[, ADJUSTABLE := PARNAME %in% c("len","sigma")]
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

# create a tempory data table with the uncertainty for each parameter
tmp <- curSysDt[ERRTYPE=="talyspar_endep"]
tmp[,DIDX:=NULL]
tmp[,IDX:=NULL]
tmp[,EN:=NULL]
tmp[,ORIGIDX:=NULL]
tmp <- unique(tmp)
stopifnot(gpDt[PARNAME=="sigma"]$EXPID == tmp$EXPID)
# They should have the same order, but should make some reordering of tmp



# define the hyperparameter limits
setkey(gpDt, IDX)
lowerLims <- rep(NA_real_, nrow(gpDt))
upperLims <- rep(NA_real_, nrow(gpDt))
lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- 0.
lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- 3
lowerLims <- lowerLims[gpDt$ADJUSTABLE]
upperLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- tmp$UNC
upperLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- 200
upperLims <- upperLims[gpDt$ADJUSTABLE]

# initial configuration
#initPars <- lowerLims + runif(length(lowerLims)) * abs(upperLims - lowerLims)
initPars <- rep(NA_real_, nrow(gpDt))
initPars[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- 0.5*tmp$UNC
initPars[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- 50
initPars <- initPars[gpDt$ADJUSTABLE]
#initPars <- upperLims

# optimize hyperparameters
# this may take a few minutes
#optRes <- optim(initPars, optfuns$logLike, optfuns$gradLogLike, method = "L-BFGS-B",
#                lower = lowerLims, upper = upperLims, control = list(fnscale = -1))

library(optimParallel)
# Setup of multicore optimization using optimparalell
nCores <- detectCores(all.tests = FALSE, logical = TRUE)
cl <- makeCluster(min(nCores,6))
setDefaultCluster(cl=cl)
#optimParallel
optRes <- optimParallel(par = initPars, 
                        fn = optfuns$logLike, 
                        gr = optfuns$gradLogLike, 
                        method = "L-BFGS-B",
                        lower = lowerLims, 
                        upper = upperLims, 
                        control = list(fnscale = -1)
)

newDts <- optfuns$getModifiedDts(optRes$par)

optExpDt <- newDts$expDt
optSysDt <- newDts$sysDt
optGpDt <- newDts$gpDt
optGpDt <- optGpDt[order(IDX)]
optGpDt[ADJUSTABLE==TRUE,INITVAL:=initPars]

# It appears that after MLO the hyper-parameter sigma get one of 
# two values:
# (1) the maximum allowed value := prior parameter uncertainty
# (2) zero!
# 

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)

print("---- optimised endep hyperparameters -----")
print(optGpDt[ADJUSTABLE==TRUE],topn=nrow(optGpDt[ADJUSTABLE==TRUE]))

