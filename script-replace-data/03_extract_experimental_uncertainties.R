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

if(length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
} else if (length(args) > 1) {
  stop("Script only accepts one argument.", call.=FALSE)
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
refParamDt <- read_object(2, "refParamDt")

#################################################
#       START OF SCRIPT
##################################################


print("-----------------------------------------------------")
print("----------------------script 03----------------------")
print("-----------------------------------------------------")

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

# SUBENTRY 22433007:
# Does not comply with the expectation of Georgs exforUncertainty package:
# the subentry reports a mixture of statisitical and systematic uncertainties in
# ERR-1: the statistical unc. varies 0.2-4%, while there are other components that
# are systematic (constant fractional uncs.). The expectation is that ERR-1 is a 
# systematic component (constant fractional/absolute unc.). The package checks that
# the numbers are within 0.1%, if not it throws an error.
# This should be fixed: Could check if ERR-S column exist (which is does not for 
# this subentry).
# (1) If it does, you "punish" the data set by adding the maximum of the 
# listed uncertainties as a systematic component.
# (2) If it does not you assume that the column is a statistical (uncorrelated) unc.

# for now: just do this manually
idx_22433007 <- which(sapply(subents, function(x) { x$ID=="22433007" }))
stopifnot(subents[[idx_22433007]]$ID == "22433007") # just double check that we change the correct subent
# change the specification of ERR-1 to ERR-S (i.e. systematic to statistical)
colnames(subents[[idx_22433007]]$DATA$TABLE) <- c("EN-ERR", "ERR-2", "EN", "DATA", "ERR-S", "ERR-T")
# make the same change in the data description, and data err
subents[[idx_22433007]]$DATA$DESCR <- c("EN-ERR", "ERR-2", "EN", "DATA", "ERR-S", "ERR-T")
names(subents[[idx_22433007]]$DATA$UNIT) <- c("EN-ERR", "ERR-2", "EN", "DATA", "ERR-S", "ERR-T")


rawUncDt <- getAvailableUncertainties(expDt, subents, dataref.col = "DATAREF")
expUncDt <- compactifyUncDt(rawUncDt)
expUncStruc <- splitUncDt(expUncDt)

sysUncDt <- expUncStruc$sysUncDt
# violation of R design philosophy: we change expDt by reference!
addStatUncToExpDt(expDt, expUncStruc$statUncDt)

# Sometimes a subentry can be incomplete in the sense that
# stat. err is not given for all points, maybe this should be removed
# but it causes trouble with the covariance matrix. For now I will just
# assume the stat. err is the maximum of 10% or the largest error present
if(any(is.na(expDt[,UNC]))) {
  cat("======================== WARNING ========================\n")
  cat("= missing stat. unc in entries:                         =\n")
  cat("=========================================================\n")
  for(curExpId in unique(expDt[is.na(UNC),EXPID])) {
    curExpDt <- expDt[EXPID==curExpId]
    print(curExpDt)
    maxUNCrel <- max(c(0.1,curExpDt[,UNC/DATA]),na.rm=TRUE)
    expDt[EXPID==curExpId & is.na(UNC),UNC:=maxUNCrel*DATA]
    cat("======================== CORRECTION ========================\n")
    print(expDt[EXPID==curExpId])
    cat("=========================================================\n")
  }
}

# safeguard against very small and very large variations
# both may cause trouble in the GLS procedure
expDt[, UNC := pmin(1000, pmax(UNC, 1))]


# model predictions mapped to dense expgrid
energyGrid <- unlist(refParamDt[PARNAME=="energy",PARVAL])
modList <- expDt[, list(SUBENT=list(createSubentStub(.BY[["REAC"]], energyGrid))), by="REAC"]
modDt <- exforHandler$extractData(modList$SUBENT, ret.values = FALSE)
modDt <- exforHandler$map(modDt, extNeedsDt, modList$SUBENT)
Smodexp <- exforHandler$getJac(modDt, extNeedsDt, modList$SUBENT)

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)
