#
# DESCRIPTION OF STEP
#
# Use the 2nd order Taylor approximation of
# the posterior to create a sample of 
# TALYS parameter sets. Perform the respective
# calculations and save the results.
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

scriptnr <- 12L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

extNeedsDt <- read_object(2, "extNeedsDt")
optParamDt <- read_object(10, "optParamDt")
needsDt <- read_object(1, "needsDt")
Sexp <- read_object(10, "Sexp")
mask <- read_object(10, "mask")
optSysDt_allpars <- read_object(10, "optSysDt_allpars")
finalPars <- read_object(11, "finalPars")
finalParCovmat <- read_object(11, "finalParCovmat")


##################################################
#       START OF SCRIPT
##################################################

# define objects to be returned
outputObjectNames <- c("allParsets", "allResults")
check_output_objects(scriptnr, outputObjectNames)

# create data table for extrapolation
extNeedsDt <- needsDt[,{
  stopifnot(all(L2 == 0) & all(L3 == 0))
  list(L1 = defineEnergyGrid(L1, energyGridrandomFiles, enPolicy="compgrid"),
       L2 = 0, L3 = 0, V1 = 0)
}, by=c("PROJECTILE", "ELEMENT", "MASS", "REAC")]

extNeedsDt[, IDX := seq_len(.N)]

# now consider also the insensitive parameters
# available for variations
optParamDt[PARNAME %in% optSysDt_allpars$PARNAME, ADJUSTABLE:=TRUE]

# set ener gygrid for random-files
parval <- as.list(optParamDt[, PARVAL])
parval[[1]] <- as.vector(energyGridrandomFiles)
optParamDt[, PARVAL:= parval]

# see step 07_tune_talyspars.R for more explanation
# about setting up the talys handler
talysHnds <- createTalysHandlers()
talys <- talysHnds$talysOptHnd
talys$setPars(optParamDt)
talys$setParTrafo(paramTrafo$fun, paramTrafo$jac)
talys$setNeeds(extNeedsDt)
talys$setSexp(Sexp)
talys$setMask(mask)
talys$setEps(0.01)

# set the seed for the random number generator
# to have a reproducible creation of TALYS parameter sets
set.seed(talysFilesSeed)

# create a sample of parameter sets
variedParsets <- sample_mvn(numTalysFiles, finalPars, finalParCovmat)

allParsets <- cbind(finalPars, variedParsets)

# perform calculations and save the result
talysHnds$remHnd$ssh$execBash(paste0("mkdir -p '", pathTalys, "'; echo endofcommand"))

allResults <- talys$fun(allParsets, applySexp = FALSE, ret.dt=FALSE, saveDir = savePathTalys)

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)
