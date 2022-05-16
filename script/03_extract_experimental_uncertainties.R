#
# DESCRIPTION OF STEP
#
# Retrieve the uncertainties from the EXFOR 
# entries and add extra uncertainties if
# something is fishy, e.g., individual
# uncertainty components do not add up to the
# total uncertainty or systematic uncertainties
# are missing. The reference calculation of the
# step 02 is used as the basis to convert
# relative uncertainties to absolute ones
# in order to avoid Peelle's pertinent puzzle.
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

scriptnr <- 3L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

subents <- read_object(1, "subents")
expDt <- read_object(1, "expDt")
extNeedsDt <- read_object(2, "extNeedsDt")

#################################################
#       START OF SCRIPT
##################################################

# define the objects that should be returned
outputObjectNames <- c("expDt", "sysUncDt",
                       "modList", "modDt", "Smodexp")
check_output_objects(scriptnr, outputObjectNames)

# augment expDt with the model prediction
Sexp <- exforHandler$getJac(expDt, extNeedsDt, subents)
setkey(expDt, IDX)
setkey(extNeedsDt, IDX)
expDt[, DATAREF := as.vector(Sexp %*% extNeedsDt$V1)]
expDt[, REFDATA := 0]

# difference between DATAREF and REFDATA
# DATAREF is used as a baseline cross section
#         to convert absolute cross sections uncertainties into
#         relative ones (we do not take the
#         experimental xs to avoid Peelle's pertinent puzzle)
# REFDATA Not used in this file. REFDATA provides the
#         predicted values associated with a reference parameter
#         set. The use of REFDATA becomes necessary for 
#         inhomogenous and/or non-linear functions that should
#         approximated by a first order Taylor polynomial

rawUncDt <- getAvailableUncertainties(expDt, subents, dataref.col = "DATAREF")
expUncDt <- compactifyUncDt(rawUncDt)
expUncStruc <- splitUncDt(expUncDt)

sysUncDt <- expUncStruc$sysUncDt
# violation of R design philosophy: we change expDt by reference!
addStatUncToExpDt(expDt, expUncStruc$statUncDt)

# safeguard against very small and very large variations
# both may cause trouble in the GLS procedure
expDt[, UNC := pmin(1000, pmax(UNC, 1))]


# model predictions mapped to dense expgrid

modList <- expDt[, list(SUBENT=list(createSubentStub(.BY[["REAC"]], energyGrid))), by="REAC"]
modDt <- exforHandler$extractData(modList$SUBENT, ret.values = FALSE)
modDt <- exforHandler$map(modDt, extNeedsDt, modList$SUBENT)
Smodexp <- exforHandler$getJac(modDt, extNeedsDt, modList$SUBENT)

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)
