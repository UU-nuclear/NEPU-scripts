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
  #source("./config/config.R")
  #stop("No config file supplied, using default file config.R", call.=FALSE)
} else if (length(args) > 1) {
  stop("Script only accepts one argument.", call.=FALSE)
} else {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}

talysHnds <- createTalysHandlers()

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



print("-----------------------------------------------------")
print("----------------------script 05----------------------")
print("-----------------------------------------------------")

# define objects to be returned
outputObjectNames <- c("fullSensDt")
check_output_objects(scriptnr, outputObjectNames)

# create the parameter transformation object
adjParNames <- refParamDt[ADJUSTABLE==TRUE,PARNAME]
for( i in 1:nrow(parRanges) ) {
  refParamDt[grepl(parRanges[i]$keyword,PARNAME),PARMIN:=parRanges[i]$min]
  refParamDt[grepl(parRanges[i]$keyword,PARNAME),PARMAX:=parRanges[i]$max]
}

# set the parameter transformation to be centered at the prior mean/mode
# the ranges of the parameters are the ones specified in the TALYS manual
paramTrafo <- parameterTransform(
                  x0 = unlist(refParamDt[ADJUSTABLE==TRUE,PARVAL]),
                  delta = refParamDt[ADJUSTABLE==TRUE,unlist(PARVAL) - PARMIN])

# We let the talys_wrapper take care of creating the jacobian
talys <- talysHnds$talysOptHnd

talys$setPars(refParamDt)
# the optimization will be undertaken using transformed
# TALYS parameters to make sure that parameters remain
# within given limits. paramTrafo is defined in config.R
talys$setParTrafo(paramTrafo$fun, paramTrafo$jac)
# define the required predictions from TALYS in order to
# map to the grid of the experimental data given in EXFOR
talys$setNeeds(extNeedsDt)
# the finite difference used to calculate numerically the
# Jacobian matrix
if(!exists("talys_finite_diff")) talys_finite_diff <- 0.01
talys$setEps(talys_finite_diff)

# We could speed up the calculation a bit here by adding a default mask:
# We know that energy dependent parameters are never sensitive to the full
# energy range, for example the parameter value at 10 MeV only affects
# the cross section at 10 MeV and above, never below. So a default mask
# could be made that limits the energy points needed for each energy
# dependent parameter. This would cut the number of calculations for these
# paramters in half.

# now we perform the calculation of the Jacobian (aka sensitivity matrix).
# Note concerning parameter transformation: values in refParamDt are assumed to
# be in the original parameter space. 
# we have two parameter spaces:
# internal: as seen by the LM algorithm, these parameters are not limited 
#           and may take on any real value in +-inf
# external: as seen by talys, these parameters are limited by the transformation
#           to the interval p0 +- delta, where p0 is the prior expectation.
# The talys_wrapper transforms the values given to it into the external parameter
# space before doing any talys calculation. ( Note that the transformation is
# constructed so that p_int = p_ext when the parameters are at their prior
# expectation values p0.) 
# talys$jac() with returnDt=TRUE will return a data table with the
# results in the external parameter space. In order to get it in the internal
# parameter space it needs to be mulitplied with paramTrafo$jac(p_int)
cat("Started calculations at", as.character(Sys.time()), "\n")
fullSensDt <- talys$jac(unlist(refParamDt[ADJUSTABLE==TRUE,PARVAL]),returnDt=TRUE)
cat("Finished calculations at", as.character(Sys.time()), "\n")

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)

