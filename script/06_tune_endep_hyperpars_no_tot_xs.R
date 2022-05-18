#
# DESCRIPTION OF STEP
#
# tune the hyperparameters of the Gaussian
# processes attached to the energy-dependent
# TALYS parameters using maximunm likelihood
# optimization.
#
# Here: exclude the total cross section to try
# to diminish the bias from the URR

#################################################
#       SCRIPT Setup
##################################################

args = commandArgs(trailingOnly=TRUE)


if(length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
} else if(length(args)>1) {
  stop("Script only accepts one argument.", call.=FALSE)
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

# remove all data on (n,tot) for the GP hyper-parameter optimization
expDt_part <- expDt[REAC!="(26-FE-56(N,TOT),,SIG)"]
expDt_part[, IDX := seq_len(.N)]

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
    gpHandler$addGP(curParname, sigma=0.1, len=2, nugget=1e-04)
}

# create a mapping matrix from the model output
# to the grid of the experimental data
Sexp <- exforHandler$getJac(expDt_part, extNeedsDt, subents) 

# prepare the experimental data
curExpDt <- copy(expDt_part) 
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
# what happens if I allow the "noise" parameter be free? => nothing: the nugget parameter gets the init value
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

# this is strange, why does the optimization of the GP hyperparameters change the experimental data? does it really?
#optExpDt <- newDts$expDt
optSysDt <- newDts$sysDt
optGpDt <- newDts$gpDt

# set the output experimental data and systematic components to the input of the same
# make the same modifications to the structure as made for the subset of data without
# the (n,tot) cross-section. This is just for the next step to get what it expects
optExpDt <- copy(expDt)
setkey(optExpDt, IDX)
# create a mapping matrix from the model output
# to the grid of the experimental data
Sexp <- exforHandler$getJac(optExpDt, extNeedsDt, subents) 
# REFDATA contains the cross sections of the reference model
# calculation. It serves as the expansion point in the
# Taylor approximation of the model entering the Bayesian
# update
optExpDt[, REFDATA := as.vector(Sexp %*% extNeedsDt$V1)]
optExpDt[, ADJUSTABLE := FALSE]


# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)

