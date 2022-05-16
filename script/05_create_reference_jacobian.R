#
# DESCRIPTION OF STEP
#
# Calculate the sensitivity matrix of TALYS
# taking all adjustable model parameters into
# account. This matrix is required for step 06
# to tune the hyperparameters of the Gaussian
# processes defined on energy-dependent TALYS
# parameters.
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

scriptnr <- 5L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

refInpList <- read_object(2, "refInpList")
refParamDt <- read_object(2, "refParamDt")
extNeedsDt <- read_object(2, "extNeedsDt") 

##################################################
#       START OF SCRIPT
##################################################

# define objects to be returned
outputObjectNames <- c("fullSensDt")
check_output_objects(scriptnr, outputObjectNames)

# sanity check: are all parameters within the restricted
# range allowed by the transformation
# (check near end in config.R for details about the transformation)
trafoPars <- paramTrafo$invfun(refParamDt[ADJUSTABLE == TRUE, unlist(PARVAL)])
origPars <- paramTrafo$fun(trafoPars)
stopifnot(all.equal(refParamDt[ADJUSTABLE==TRUE, unlist(PARVAL)], origPars))

# create the inputs for the calculation of the Jacobian (aka sensitivity matrix)
# Note concerning parameter transformation: values in refParamDt are assumed to
# be in the original parameter space. Also generated input lists in jacInputsDt
# contain untransformed parameters. However, the eps specification refers to 
# adjustments in the transformed! parameter space
jacInputsDt <- createInputsForJacobian(refParamDt, extNeedsDt, eps = 0.01, trafo = paramTrafo)

###############################
# perform the calculation
##############################

talysHnds <- createTalysHandlers()
talysHnd <- talysHnds$talysHnd

runObj <- talysHnd$run(jacInputsDt$inputs, jacInputsDt$outspecs)
# save the information about the job for 
# later recovery if something goes wrong
save_output_objects(scriptnr, "runObj", overwrite)

cat("Started calculations at", as.character(Sys.time()), "\n")
cat("Waiting for termination...\n")
while (talysHnd$isRunning(runObj)) { 
    Sys.sleep(pollTime)
}
jacRes <- talysHnd$result(runObj)
talysHnds$clustHnd$closeCon()

jacInputsDt$outspecs <- lapply(jacRes, function(x) x$result)
fullSensDt <- computeJacobian(jacInputsDt, drop0=TRUE)

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)

