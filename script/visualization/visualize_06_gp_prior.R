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

library(ggplot2)
library(stringr)
library(data.table)
library(Matrix)

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

outdataPathRun <- outdataPath
plotPath <- paste0(outdataPathRun, "/plots")


subents <- read_object(1, "subents", outdata_path = outdataPathRun)
modList <- read_object(3, "modList", outdata_path = outdataPathRun)
optSysDt <- read_object(7, "optSysDt", outdata_path = outdataPathRun)
optGpDt <- read_object(7, "optGpDt", outdata_path = outdataPathRun)
optExpDt <- read_object(7, "optExpDt", outdata_path = outdataPathRun)
mapAssignment <- read_object(7, "mapAssignment", outdata_path = outdataPathRun)
sysCompHandler <- read_object(7, "sysCompHandler", outdata_path = outdataPathRun)
origSysDt <- read_object(4, "origSysDt", outdata_path = outdataPathRun)
updSysDt <- read_object(4, "updSysDt", outdata_path = outdataPathRun)
expDt <- read_object(3, "expDt", outdata_path = outdataPathRun)
yexp <- read_object(7, "yexp", outdata_path = outdataPathRun)
optRes <- read_object(7, "optRes", outdata_path=outdataPathRun)
normHandler <- createSysCompNormHandler("DATAREF")
normHandler$addSysUnc("EXPID", "", 0, 0, TRUE)
fitResult <- as.numeric(as.vector(optRes$fn))

sysCompHandler2 <- createSysCompHandler()
sysCompHandler2$addHandler(normHandler)

S <- sysCompHandler2$map(expDt, origSysDt, ret.mat = TRUE)
origX <- sysCompHandler2$cov(origSysDt, ret.mat = TRUE)
updX <- sysCompHandler2$cov(updSysDt, ret.mat = TRUE)
statUnc <- getDt_UNC(expDt)

origUnc <- sqrt(statUnc^2 + diag(S %*% origX %*% t(S))) 
updUnc <- sqrt(statUnc^2 + diag(S %*% updX %*% t(S)))


optExpDt[, LMFIT := fitResult]
optExpDt[, UPDUNC := updUnc]

optExpDt <- optExpDt[order(-L1),.SD,by=REAC]

#optSysDtGpObs <- optSysDt[grepl("REACEXP-", EXPID) & EXPID == "REACEXP-N,TOT-01", ][, IDX := seq_len(.N)]
#Sk <- sysCompHandler$map(optExpDt[REAC=="(26-FE-56(N,TOT),,SIG)"][, IDX := seq_len(.N)], optSysDtGpObs, ret.mat = TRUE)
optSysDtGpObs <- optSysDt[grepl("REACEXP-", EXPID)][, IDX := seq_len(.N)]
Sk <- sysCompHandler$map(optExpDt[, IDX := seq_len(.N)], optSysDtGpObs, ret.mat = TRUE)

K <- sysCompHandler$cov(optSysDtGpObs, optGpDt[grepl("REACEXP-", EXPID), ][, IDX := seq_len(.N)], ret.mat = TRUE)
statUnc <- getDt_UNC(optExpDt)

gpObsPrior <- sqrt(diag(Sk %*% K %*% t(Sk)))


optExpDt[, GPPRIOR := gpObsPrior]
#expDtBay[, ORIGUNC := optExpDt$UNC]
#expDtBay[, ORIGDATA := optExpDt$DATA]

#optExpDtShort <- optExpDt[REAC=="(26-FE-56(N,2N)26-FE-55,,SIG)"]



plot_gp_prior(optExpDt, optGpDt)
