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

# This script is a try to change the GP covariance function by altering 
# the sysDt data.table. The energies of the steering points are specified
# in paramDt these are the points that talys uses to do the calculation
# Normally the energies in sysDt are the same, but if we set them 
# to something else it should alter how the covariance between the points
# are calculated.

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
refParamDt <- read_object(2, "refParamDt")
extNeedsDt <- read_object(2, "extNeedsDt")
modList <- read_object(3, "modList")
expDt <- read_object(3, "expDt")
updSysDt <- read_object(4, "updSysDt")
fullSensDt <- read_object(5, "fullSensDt")
Sexp <- read_object(5,"Sexp")
optParamDt <- read_object(5, "optParamDt")

source("script-Cr/MLOwithPrior.R")
library(parallel)
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

# Here is where the 'magic' happens
# I replace the EN values by their index
energies <- unique(curSysDt[ERRTYPE=="talyspar_endep",EN])
curSysDt[,EN:=match(EN,energies),by=EN]

# set up Gaussian processes for energy dependent parameters
gpHandler <- createSysCompGPHandler()
parnames_endep <- unique(curSysDt[ERRTYPE == "talyspar_endep", EXPID])
for (curParname in parnames_endep) {
    gpHandler$addGP(curParname, curSysDt[EXPID==curParname,UNC][1], 100, 1e-4)
}

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

# set hyper-parameters for energy dependent parameters that are not adjustable
# to not be part of the hyper-parameter optiization
adjustable_endep_par_names <- optParamDt[ADJUSTABLE==TRUE]$PARNAME[grepl("\\(.+\\)",optParamDt[ADJUSTABLE==TRUE]$PARNAME)]
adjustable_endep_par_names <- unique(str_remove(adjustable_endep_par_names,"\\(.+\\)"))
gpDt[!(str_remove(EXPID,"TALYS-") %in% adjustable_endep_par_names), ADJUSTABLE:=FALSE]

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
gpDt[PARNAME=="sigma",PARUNC:=tmp$UNC]
# They should have the same order, but should make some reordering of tmp



# define the hyperparameter limits
setkey(gpDt, IDX)
lowerLims <- rep(NA_real_, nrow(gpDt))
upperLims <- rep(NA_real_, nrow(gpDt))
lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- 0.1*gpDt[ADJUSTABLE==TRUE & PARNAME=="sigma",PARUNC]
lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- 6
upperLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- gpDt[ADJUSTABLE==TRUE & PARNAME=="sigma",PARUNC]
upperLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- 3*length(energyGridForParams)

# prior on the hyper-parameter
priorExpectation <- rep(NA_real_, nrow(gpDt))
priorExpectation[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- gpDt[ADJUSTABLE==TRUE & PARNAME=="sigma",PARUNC]
priorExpectation[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- length(energyGridForParams)

priorUnc <- rep(NA_real_, nrow(gpDt))
priorUnc[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- 0.5*gpDt[ADJUSTABLE==TRUE & PARNAME=="sigma",PARUNC]
priorUnc[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- length(energyGridForParams)

#########################

lowerLims <- lowerLims[gpDt$ADJUSTABLE]
upperLims <- upperLims[gpDt$ADJUSTABLE]
priorExpectation <- priorExpectation[gpDt$ADJUSTABLE]
priorUnc <- priorUnc[gpDt$ADJUSTABLE]
initPars <- priorExpectation

# we assume no prior correlations
priorCovMat <- Diagonal(length(priorUnc),x=priorUnc^2)

# set the prior
optfuns$setPrior(priorExpectation,priorCovMat)


library(optimParallel)
# Setup of multicore optimization using optimparalell
nCores <- detectCores(all.tests = FALSE, logical = TRUE)
cl <- makeCluster(min(nCores,20))

cat("number of cores used for optimParallel ", min(nCores,20))

setDefaultCluster(cl=cl)
dummy <- clusterEvalQ(cl, c(library(data.table)))
clusterExport(cl, c("optfuns","logLike","gradLogLike"), 
              envir=environment())
#optimParallel
cat("Started optimization at", as.character(Sys.time()), "\n") 
optRes <- optimParallel(par = initPars, 
                        fn = optfuns$logPost, 
                        gr = optfuns$gradLogpost, 
                        method = "L-BFGS-B",
                        lower = lowerLims, 
                        upper = upperLims, 
                        control = list(fnscale = -1,maxit=300)
)
cat("Finished optimization at", as.character(Sys.time()), "\n")

newDts <- optfuns$getModifiedDts(optRes$par)

optExpDt <- newDts$expDt
optSysDt <- newDts$sysDt
optGpDt <- newDts$gpDt
optGpDt <- optGpDt[order(IDX)]
optGpDt[ADJUSTABLE==TRUE,INITVAL:=initPars]

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)

print("---- optimised endep hyperparameters -----")
print(optGpDt[ADJUSTABLE==TRUE],topn=nrow(optGpDt[ADJUSTABLE==TRUE]))
