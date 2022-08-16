#################################################
#       SCRIPT Setup
##################################################

args = commandArgs(trailingOnly=TRUE)

if (length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}

scriptnr <- 7L
overwrite <- TRUE

outputObjectNames <- c("optExpDt", "optSysDt", "optGpDt", "mapAssignment", "reacHandlerGPobs", "sysCompHandler", 
                       "Xk", "S0k")

# Retrieve the information on systematic error components from step06 and add 
# the GPobs component. Then, using MLO find the best hyperparameters

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

optRes <- read_object(7,"optRes")
sysCompHandler <- read_object(6,"sysCompHandler")
optExpDt <- read_object(6,"optExpDt")
optSysDt <- read_object(6,"optSysDt")
optGpDt <- read_object(6,"optGpDt")
modDt <- read_object(3, "modDt") # this is the prior talys model mapped to the energies of the experiments
subents <- read_object(1, "subents")
modList <- read_object(3, "modList")
expDt <- read_object(3, "expDt")

# Set up the prior on observable GP

reacHandlerGPobs <- createSysCompReacHandler(c(subents, modList$SUBENT))
reacHandlerGPobs$addMap("pw", pwMap)
gpObsHandler <- createSysCompGPHandler()

gpEnGridLength <- 200
reacs <- expDt[ , .N, by=REAC][N>1 ,REAC]
for(curReac in reacs){
  if(expDt[REAC==curReac, .N] <= gpEnGridLength){
    curEnGrid <- sort(expDt[REAC==curReac, L1])
  } else if(length(seq(getThresEn(curReac, modDt, defaultThresEn), maxExpEn, by = 0.1)) > gpEnGridLength){
    curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), maxExpEn, length.out = gpEnGridLength)
  } else{
    curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), maxExpEn, by = 0.1)
  }
  
  curUncs <- c(0, 0, rep(2, length(curEnGrid)-2))
  reacHandlerGPobs$assignMapToReac("pw", curReac,
                                   vals = rep(0, length(curEnGrid)),
                                   uncs =  rep(1e4, length(curEnGrid)),
                                   opts = list(ens = curEnGrid,
                                               order = 1, outsideZero = TRUE))
  
}

curSysDtGPobs <- reacHandlerGPobs$createSysDtGpObs()

mapAssignment <- reacHandlerGPobs$getMapAssignment()[,data.table(REAC, EXPID)]

exactGPObsSysDt <- copy(expDt) 
setnames(exactGPObsSysDt, "L1", "EN")
setnames(exactGPObsSysDt, "EXPID", "EXPIDOLD")
exactGPObsSysDt <- exactGPObsSysDt[mapAssignment,on="REAC"]
exactGPObsSysDt[,GPTYPE:="sqrexp"]

reactions <- exactGPObsSysDt[ , .N, by=EXPID][N>1 ,EXPID]

for(curReac in reactions){
  gpObsHandler$addGP(expid=curReac, sigma=1, len=0.1, nugget=1e-3)
}

gpObsDt <- gpObsHandler$createGPDt()

gpObsDt[,ADJUSTABLE := FALSE]
gpObsDt[EXPID %in% exactGPObsSysDt[,.N, by=EXPID][N>1,EXPID],  ADJUSTABLE := PARNAME %in% c("sigma","len")]
#gpObsDt[EXPID %in% exactGPObsSysDt[,.N, by=EXPID][N>1,EXPID],  ADJUSTABLE := PARNAME %in% c("sigma","len","nugget")]
gpObsDt[EXPID %in% exactGPObsSysDt[,.N, by=EXPID][N==1,EXPID] & PARNAME %in% c("sigma","len"),  PARVAL := 1]
gpObsDt[EXPID %in% exactGPObsSysDt[,.N, by=EXPID][N==1,EXPID] & PARNAME %in% c("nugget"),  PARVAL := 1]
gpObsDt[, IDX := seq_len(.N)]

gpObsHandler$updateSysDt(curSysDtGPobs)
curSysDtGPobs[,IDX:=seq(1,.N)]
curSysDtGPobs[,REFDATA:=rep(0,.N)]
curSysDtGPobs[,ADJUSTABLE:=rep(FALSE,.N)]

# copy the systematics data table from step06 and add the information just created on the GPobs
curSysDt <- copy(optSysDt)
curSysDt <- rbind(curSysDt, curSysDtGPobs, fill=TRUE)
curSysDt[, IDX := seq_len(.N)]

# copy the data table with gp-hyperparameters from step 06 and add the new ones for Gp in observable
gpDt <- copy(optGpDt)
gpDt$ADJUSTABLE <- FALSE # don't re-optimise the hyperpars in parameter domain
gpDt <- rbind(gpDt, gpObsDt, fill=TRUE)
gpDt[, IDX := seq_len(.N)]

# register the GPobs handler to the global systematic error handler
sysCompHandler$addHandler(reacHandlerGPobs)

optExpDt[, REFDATA := optRes$fn] # set reference data to the result of the LM algorithm
# setup optimization specification
optfuns <- createMLOptimFuns()
optfuns$setDts(optExpDt, curSysDt, gpDt,
               sysCompHandler = sysCompHandler)


# define the hyperparameter limits
setkey(gpDt, IDX)
lowerLims <- rep(NA_real_, nrow(gpDt[ADJUSTABLE==TRUE]))
upperLims <- rep(NA_real_, nrow(gpDt[ADJUSTABLE==TRUE]))

lowerLims[gpDt[ADJUSTABLE==TRUE]$PARNAME=="sigma"] <- 1e-06
upperLims[gpDt[ADJUSTABLE==TRUE]$PARNAME=="sigma"] <- 2000
lowerLims[gpDt[ADJUSTABLE==TRUE]$PARNAME=="len"] <- 1e-06
upperLims[gpDt[ADJUSTABLE==TRUE]$PARNAME=="len"] <- 5

gpDt[ADJUSTABLE==TRUE, INITVAL := PARVAL]

# optimize hyperparameters
# this may take a few minutes
library(optimParallel)
# Setup of multicore optimization using optimparalell
nCores <- detectCores(all.tests = FALSE, logical = TRUE)
cl <- makeCluster(6)
setDefaultCluster(cl=cl)
#optimParallel
optRes <- optimParallel(par = gpDt[ADJUSTABLE==TRUE, INITVAL], 
                        fn = optfuns$logLike, 
                        gr = optfuns$gradLogLike, 
                        method = "L-BFGS-B",
                        lower = lowerLims, 
                        upper = upperLims, 
                        control = list(fnscale = -1)
)

newDts <- optfuns$getModifiedDts(optRes$par)

optExpDt <- newDts$expDt
optSysDt <- newDts$sysDt
optGpDt <- newDts$gpDt

# construct the matrix to map from systematic error components
# of the experiments to the measurement points
S <- sysCompHandler$map(optExpDt, optSysDt, ret.mat = TRUE)
# construct the covariance matrix containing both the 
# covariances of the systematic experimental errors and
# the covariances for the model parameters. This matrix is diagonal
# for the portions relating to experiments and energy-independent
# parameters. The blocks associated with energy-dependent TALYS
# parameters are constructed by a Gaussian process and contain
# therefore correlations.
P <- sysCompHandler$cov(optSysDt, optGpDt, ret.mat = TRUE)
# The Levenberg-Marquardt routine assumes that systematic
# components in optSysDt are only related to experiments
# and not to model parameters. The latter were introduced
# in step 06 to optimize the hyperparameters of the GPs
# associated with energy-dependent TALYS parameters.
# Therefore we need to remove them now.
setkey(optSysDt, IDX)
setkey(optExpDt, IDX)
expSel <- optSysDt[, !grepl("TALYS-", EXPID)]
talysSel <- !expSel


Xk <- P[expSel, expSel] 
S0k <- S[, expSel]
# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)

print("-------------------------------------------------")
print("- result of MLO optimization of hyperparameters -")
print(optGpDt[grepl("REACEXP",EXPID)])