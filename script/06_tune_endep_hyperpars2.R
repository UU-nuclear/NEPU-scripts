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
min_sigma <- 0.1
max_sigma <- 0.5
min_length <- 3
max_length <- 30
setkey(gpDt, IDX)
lowerLims <- rep(NA_real_, nrow(gpDt))
upperLims <- rep(NA_real_, nrow(gpDt))
lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- min_sigma
lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- min_length
lowerLims <- lowerLims[gpDt$ADJUSTABLE]
upperLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- max_sigma
upperLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- max_length
upperLims <- upperLims[gpDt$ADJUSTABLE]

# initial configuration
#initPars <- lowerLims + runif(length(lowerLims)) * abs(upperLims - lowerLims)
#initPars <- lowerLims
#initPars <- upperLims
initPars <- rep(NA_real_, nrow(gpDt))
initPars[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- 0.5
initPars[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- 5
initPars <- initPars[gpDt$ADJUSTABLE]

# optimize hyperparameters
# this may take a few minutes
optRes <- optim(initPars, optfuns$logLike, optfuns$gradLogLike, method = "L-BFGS-B",
                lower = lowerLims, upper = upperLims, control = list(fnscale = -1))
# this actually maximizes since fnscale = -1, so i should find the maximum in the sampling
#
#newDts <- optfuns$getModifiedDts(optRes$par)
#
#optExpDt <- newDts$expDt
#optSysDt <- newDts$sysDt
#optGpDt <- newDts$gpDt
#
## save the needed files for reference
#save_output_objects(scriptnr, outputObjectNames, overwrite)
n_sigmas <- length(lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"])
n_lengths <- length(lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"])

n_adjustable_hyper_pars <- n_sigmas + n_lengths
hyper_par_vector <- rep(0, n_adjustable_hyper_pars)
sum_of_weights <- 0

# Because the logLikelihood is very small I need to normalize it by adding an offse
# otherwise I get 0 for all samples
LogLikelihood_init <- optfuns$logLike(initPars)

max_n_iterations <- 10
for(i in 1:max_n_iterations) {
  print(paste0("i = ",i))
  # sample the likelihood function
  n_samples <- 1000

  sigma_samples <- runif(n_sigmas*n_samples,min_sigma,max_sigma)
  length_samples <- runif(n_lengths*n_samples,min_length,max_length)

  hyper_par_samples <- rep(NA_real_, nrow(gpDt)*n_samples)
  hyper_par_samples[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- sigma_samples
  hyper_par_samples[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- length_samples
  hyper_par_samples <- hyper_par_samples[!is.na(hyper_par_samples)] # removing the ones that are not adjustable
  hyper_par_samples <- matrix(hyper_par_samples,ncol=n_samples)
  Likelihood_samples <- exp(apply(hyper_par_samples,2,FUN=optfuns$logLike) - LogLikelihood_init)
  print("--- Likelihood_samples ---")
  print(Likelihood_samples)

  # Likelihood_samples is vector with likelihoods each entry correspond to a column in hyper_par_samples, so for each row
  # hyper_par_samples I want the weighted average, where the weight of each column is Likelihood_samples
  curr_hyper_par_vector <- apply(hyper_par_samples,1,FUN=weighted.mean,w=Likelihood_samples)
  curr_sum_of_weights <- sum(Likelihood_samples)
  print("--- curr_hyper_par_vector ---")
  print(curr_hyper_par_vector)
  if(i>1) {
    avg_hyper_par_vector <- (hyper_par_vector+curr_sum_of_weights*curr_hyper_par_vector)/(sum_of_weights+curr_sum_of_weights)
    former_avg_hyper_par_vector <- hyper_par_vector/sum_of_weights
    diff <- avg_hyper_par_vector - former_avg_hyper_par_vector

    print("--- avg_hyper_par_vector ---")
    print(avg_hyper_par_vector)

    print("--- former_avg_hyper_par_vector ---")
    print(former_avg_hyper_par_vector)

    print("--- rel. diff ---")
    print(diff/avg_hyper_par_vector)

    if(all(abs(diff/avg_hyper_par_vector) < 1.e-03)) {
      # this convergence condition is likely to be fulfilled as soon as 
      # the samples happen to fall at small likelyhoods, I need to make sure that the 
      # current sample is representative!
      sum_of_weights <- sum_of_weights + curr_sum_of_weights
      hyper_par_vector <- hyper_par_vector + curr_sum_of_weights*curr_hyper_par_vector
      print(paste("converged after",i,"iterations"))
      break
    }
  }
  sum_of_weights <- sum_of_weights + curr_sum_of_weights
  hyper_par_vector <- hyper_par_vector + curr_sum_of_weights*curr_hyper_par_vector
  if(i==max_n_iterations) {
    print("the sampling of hyperparameters didn't converge")
  }
}

hyper_par_vector <- hyper_par_vector/sum_of_weights