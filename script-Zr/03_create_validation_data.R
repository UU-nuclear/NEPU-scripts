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
#check_output_objects(scriptnr, outputObjectNames)


outdataPath_original <- outdataPath

# augment expDt with the model prediction
Sexp <- exforHandler$getJac(expDt, extNeedsDt, subents)
setkey(expDt, IDX)
setkey(extNeedsDt, IDX)
expDt[, DATAREF := as.vector(Sexp %*% extNeedsDt$V1)]
expDt[, REFDATA := 0]

# divide the full data set into 5 sets
# expDt[, DATASET := sample.int(5,.N,replace=TRUE)]

# split the data randomly into 5 equally sized data.tables 
expDt[,DATASET:=sample(x=rep(1:5,times=ceiling(.N/5)),size=.N,replace=FALSE)]

expDt_orig <- copy(expDt)

for(DATASET_IDX in 1:5) {
  expDt <- expDt_orig[DATASET!=DATASET_IDX]
  expDt[,IDX:=seq_len(.N)]

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

  outdataPath <- file.path(outdataPath_original, paste0("validation",DATASET_IDX))

  dir.create(outdataPath, recursive=TRUE, showWarnings=FALSE)
  check_output_objects(scriptnr, outputObjectNames)
  save_output_objects(scriptnr, outputObjectNames, overwrite)
}
