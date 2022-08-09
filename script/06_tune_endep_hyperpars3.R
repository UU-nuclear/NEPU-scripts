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

library(parallel)

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
min_sigma <- 0.1
max_sigma <- 0.5
min_length <- 3
max_length <- 10
setkey(gpDt, IDX)
lowerLims <- rep(NA_real_, nrow(gpDt))
upperLims <- rep(NA_real_, nrow(gpDt))
lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- min_sigma
lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- min_length
lowerLims <- lowerLims[gpDt$ADJUSTABLE]
upperLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- max_sigma
upperLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- max_length
upperLims <- upperLims[gpDt$ADJUSTABLE]


# optimize hyperparameters
# randomly select different starting points from a uniform prior, optimize for each starting point
# and select the one with the maximum likelihood

n_sigmas <- length(lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"])
n_lengths <- length(lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"])

# sample hyperparameter initial values
n_samples <- 64

sigma_samples <- runif(n_sigmas*n_samples,min_sigma,max_sigma)
length_samples <- runif(n_lengths*n_samples,min_length,max_length)

hyper_par_samples <- rep(NA_real_, nrow(gpDt)*n_samples)
hyper_par_samples[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- sigma_samples
hyper_par_samples[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- length_samples
hyper_par_samples <- hyper_par_samples[!is.na(hyper_par_samples)] # removing the ones that are not adjustable

list_of_samples <- as.list(split(hyper_par_samples,  cut(seq_along(hyper_par_samples), n_samples, labels = FALSE)))

num_cores <- min(detectCores(),n_samples)
optRes_list <- mclapply(list_of_samples, FUN=optim, fn = optfuns$logLike, gr = optfuns$gradLogLike,
                          method = "L-BFGS-B",lower = lowerLims, upper = upperLims, control = list(fnscale = -1),
                          mc.cores = num_cores)

optRes <- optRes_list[[1]]
for(res in optRes_list) {
  if(res$value > optRes$value) {
    optRes <- res
  }
}

# the weighted averaged result
optPar_matrix <- matrix(NA_real_, nrow = length(optRes_list[[1]]$par), ncol = n_samples)
weight_vector <- rep(NA_real_, n_samples)
idx <- 1
for(res in optRes_list) {
  optPar_matrix[,idx] <- res$par
  weight_vector[idx] <- exp(res$value - optRes$value)

  idx <- idx + 1
}

optPar_avg <- apply(optPar_matrix,1,FUN=weighted.mean,w=weight_vector)

#newDts <- optfuns$getModifiedDts(optRes$par) # using the par-vector with the maximum marginal likelihood
newDts <- optfuns$getModifiedDts(optPar_avg) # using the average par-vector weighted with the marginal likelihood

# I would expect that hyperparameters which are constrained by data will obtain the same value,
# no matter the starting point. This seems to be shown by the data as well. Fot the hyperparameters
# that are not constrained by the data the result of the MLO lands on the initial value. Hence,
# by taking the weighted average I should get the optimum for the hyperparameters that are
# constrained, while those that are not should obtain a value close to the prior (uniform) average

optExpDt <- newDts$expDt
optSysDt <- newDts$sysDt
optGpDt <- newDts$gpDt

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)