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

##################################################################

rebuildGauss <- FALSE
if (any(file.exists("modelFullData.rds"))  & !rebuildGauss){
  modelFullData <- readRDS("modelFullData.rds")
} else {
  nvar <- 1
  
#  modelFullData <- expDtTot[, 
#                            mleHetGP(X = L1, 
#                                     Z = DATA, 
#                                     lower = rep(0.1, nvar), 
#                                     upper = rep(50, nvar),
#                                     covtype = "Gaussian",
#                                     init = list(
#                                       beta0 = expDtTot$DATAREF,
#                                       theta = 2)
#                            )
#                            ]
  modelFullData <- expDt[, 
                            mleHetGP(X = L1, 
                                     Z = DATA, 
                                     lower = rep(0.1, nvar), 
                                     upper = rep(50, nvar),
                                     covtype = "Gaussian",
                                     init = list(
                                       beta0 = expDt$DATAREF,
                                       theta = 2)
                            )
                            ]
  modelFullData <- rebuild(modelFullData, robust = TRUE)
  saveRDS(modelFullData, "modelFullData.rds")
}
##############################################################

expDtTot <- expDt[REAC=="(26-FE-56(N,TOT),,SIG)"]

n = 0
grid <- seq(min(expDtTot$L1)-0.1, max(expDtTot$L1)+0.1, length.out=ceiling(max(expDtTot$L1-min(expDtTot$L1))))
weightDt <- data.table(IDX = expDtTot$IDX)

for (i in 1:length(grid)) {
  # For each bin grid, attache a waight to the data point 
  # according to the number of data points in that bin.
  weightDt[IDX %in% expDtTot[L1 >= grid[i] & L1 < grid[i+1],  IDX ], PROBDENS:= expDtTot[L1 >= grid[i] & L1 < grid[i+1],  1/.N ] ] 
  
  # Count the bins that contain non-zero number of events
  n = n + expDtTot[L1 >= grid[i] & L1 < grid[i+1],  ifelse(.N >0, 1, 0) ]
}

# Normalize the weights accoding to the number of non-zero bins
weightDt[,PROBDENS := PROBDENS/n]
expDtSampleTmp <- expDtTot[IDX %in% expDtTot[, sample(IDX, 400, prob = weightDt$PROBDENS)]]

xgrid <- as.matrix(expDtSampleTmp$L1)
predictions <- predict(x = xgrid, object =  modelFullData)


ggr <- ggplot(NULL) + theme_bw()    
#ggr <- ggr +  ggdata(expDt[REAC=="(26-FE-56(N,TOT),,SIG)"]) 
ggr <- ggr + geom_point(data = expDtTot,   aes(x = L1, y = DATA), size = 0.1, alpha=0.5)    

#
ggr <- ggr + geom_errorbar(aes(x = xgrid, ymin = predictions$mean - sqrt(predictions$sd2+ predictions$nugs), ymax = predictions$mean + sqrt(predictions$sd2+ predictions$nugs)), color="tomato1", alpha=0.4, size=0.5)     

#ggr <- ggr + geom_errorbar(aes(x = xgrid, ymin = qnorm(0.32, predictions$mean, sqrt(predictions$sd2+ predictions$nugs)), ymax = qnorm(0.68, predictions$mean, sqrt(predictions$sd2+ predictions$nugs))), color="red", alpha=0.8, size=1)    
ggr <- ggr + geom_point(aes(x = xgrid, y = predictions$mean), color="tomato1", size = 1)  
ggr <- ggr + geom_line(aes(x = xgrid, y = predictions$mean), linetype = "dashed", color="black", size = 0.5) 
print(ggr)
plot(expDtSampleTmp$L1, expDtSampleTmp$DATA)
#expDtSampleTmp <- expDtSampleTmp[sort(L1)]
plot(expDtSampleTmp$L1, predictions$mean)
#expDtSampleTmp <- expDtSampleTmp[sort(L1)]
expDtSampleTmp[, DATA := predictions$mean]
expDtSampleTmp[, UNC := sqrt(predictions$sd2+ predictions$nugs)]

plot(expDtSampleTmp$L1, expDtSampleTmp$DATA)

expDtTmp <- rbind(expDt[REAC!="(26-FE-56(N,TOT),,SIG)"], expDtSampleTmp)
expDtTmp <- expDtTmp[sort(L1), .SD, by=REAC][, IDX := seq_len(.N)]
expDt <- expDtTmp
########################################## WOOP
# model predictions mapped to dense expgrid

modList <- expDt[, list(SUBENT=list(createSubentStub(.BY[["REAC"]], energyGrid))), by="REAC"]
modDt <- exforHandler$extractData(modList$SUBENT, ret.values = FALSE)
modDt <- exforHandler$map(modDt, extNeedsDt, modList$SUBENT)
Smodexp <- exforHandler$getJac(modDt, extNeedsDt, modList$SUBENT)

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)
