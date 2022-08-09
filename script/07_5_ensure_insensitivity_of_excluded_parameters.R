#
# DESCRIPTION OF STEP
#
# 
#

#################################################
#       SCRIPT Setup
##################################################

args = commandArgs(trailingOnly=TRUE)


if (length(args)==0) {
  #source("./config/config.R")
  #stop("No config file supplied, using default file config.R", call.=FALSE)
} else if (length(args) > 1) {
  stop("Script only accepts one argument.", call.=FALSE)
} else {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}

#################################################
#       SCRIPT PARAMETERS
##################################################

scriptnr <- 7L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

subents <- read_object(1, "subents")
refParamDt <- read_object(2, "refParamDt")
extNeedsDt <- read_object(2, "extNeedsDt")
modList <- read_object(3, "modList")
#fullSensDt <- read_object(5, "fullSensDt") 
optExpDt <- read_object(6, "optExpDt")
optSysDt <- read_object(6, "optSysDt")
optGpDt <- read_object(6, "optGpDt")

optParamDt <- read_object(7, "optParamDt")
finalParamDt <- read_object(8, "finalParamDt")

# define objects to be returned
outputObjectNames <- c("fullSensDtPost")
check_output_objects(scriptnr, outputObjectNames)

# create a data.table witht he same struct as refParamDt but with the final parameter values and uncertainties
refParamDt_post <- copy(refParamDt)
refParamDt_post$PARVAL <- finalParamDt$POSTVAL
refParamDt_post$PARUNC <- finalParamDt$POSTUNC

# create the inputs for the calculation of the Jacobian (aka sensitivity matrix)
# Note concerning parameter transformation: values in refParamDt are assumed to
# be in the original parameter space. Also generated input lists in jacInputsDt
# contain untransformed parameters. However, the eps specification refers to 
# adjustments in the transformed! parameter space
#jacInputsDt <- createInputsForJacobian(refParamDt, extNeedsDt, eps = 0.001, trafo = paramTrafo)
#jacInputsDt <- createInputsForJacobian(optParamDt, extNeedsDt, eps = 0.001, trafo = paramTrafo)
jacInputsDt <- createInputsForJacobian(refParamDt_post, extNeedsDt, eps = 0.001, trafo = paramTrafo)

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
jacRes <- talysHnd$result(runObj)
cat("Finished calculations at", as.character(Sys.time()), "\n")
#talysHnds$clustHnd$closeCon()

save_output_objects(scriptnr, "jacRes", overwrite)

jacInputsDt$outspecs <- lapply(jacRes, function(x) x$result)
fullSensDtPost <- computeJacobian(jacInputsDt, drop0=TRUE)

# save the needed files for reference
save_output_objects(scriptnr, "fullSensDtPost", overwrite)
