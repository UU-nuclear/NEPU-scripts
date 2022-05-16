#
# DESCRIPTION OF STEP
#
# tune the hyperparameters of the Gaussian
# processes attached to the energy-dependent
# TALYS parameters using maximunm likelihood
# optimization.
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

scriptnr <- 6L
overwrite <- TRUE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

subents <- read_object(1, "subents")
modList <- read_object(3, "modList")
refParamDt <- read_object(2, "refParamDt")
extNeedsDt <- read_object(2, "extNeedsDt")
expDt <- read_object(3, "expDt")
updSysDt <- read_object(4, "updSysDt")
fullSensDt <- read_object(5, "fullSensDt") 
modDt <- read_object(3, "modDt")
priorDt <- read_object(4, "priorDt")
##################################################
#       START OF SCRIPT
##################################################

# define objects to be returned
outputObjectNames <- c("optExpDt", "optSysDt", "optGpDt", "mapAssignment", "reacHandlerGPobs", "sysCompHandler")
check_output_objects(scriptnr, outputObjectNames)

# prepare the data table with systematic components
# this comprises both systematic experimental 
# errors and the model parameters

talysHandler <- createSysCompModelHandler()
talysHandler$setRef(extNeedsDt, fullSensDt, refParamDt,
                    exforHandler, c(subents, modList$SUBENT))
talysHandler$setPrior(refParamDt)
modelSysDt <- talysHandler$createSysDt()

setkey(updSysDt, IDX)
curSysDt <- rbind(updSysDt, modelSysDt, fill=TRUE)
curSysDt <- curSysDt[!grepl("REACEXP-", EXPID)]
curSysDt[, IDX := seq_len(.N)]
curSysDt[, ADJUSTABLE := FALSE]

# set up Gaussian processes for energy dependent parameters
gpHandler <- createSysCompGPHandler()
parnames_endep <- unique(curSysDt[ERRTYPE == "talyspar_endep", EXPID])
for (curParname in parnames_endep) {
    gpHandler$addGP(curParname, 0.1, 3, 1e-4)
}

# create a mapping matrix from the model output
# to the grid of the experimental data
Sexp <- exforHandler$getJac(expDt, extNeedsDt, subents) 

# prepare the experimental data
curExpDt <- copy(expDt)
setkey(curExpDt, IDX)
setkey(extNeedsDt, IDX)
# REFDATA contains the cross sections of the reference model
# calculation. It serves as the expansion point in the
# Taylor approximation of the model entering the Bayesian
# update
curExpDt[, REFDATA := as.vector(Sexp %*% extNeedsDt$V1)]
curExpDt[, ADJUSTABLE := FALSE]

gpDt <- gpHandler$createGPDt()
gpHandler$updateSysDt(curSysDt)

# define which GP hyperparameters are adjustable
gpDt[, ADJUSTABLE := PARNAME %in% c("sigma","len")]
gpDt[, IDX := seq_len(.N)]

# initialize handler for systematic experimental uncertainties
# TODO: It should not be needed to define the systematic error
#       components during the setup. The category (column name)
#       should be enough.
normHandler <- createSysCompNormHandler("DATAREF")
tmpDt <- updSysDt[grepl("EXPID-", EXPID)]
expIds <- sub("EXPID-", "", tmpDt$EXPID)
normHandler$addSysUnc("EXPID", expIds, 0, tmpDt$UNC,
                      rel = (tmpDt$ERRTYPE == "sys-rel"))

# Set up the pr prior on observable

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


  #prior <- ifelse(curReac %in% priorDt$REAC,  as.numeric(priorDt[REAC==curReac, PRIOR]), 10)
  
  #curUncs <- c(1e4, 1e4, rep(prior, length(curEnGrid)-2)) # order = 3
  #reacHandlerGPobs$assignMapToReac("pw", curReac,
  #                                 vals = rep(0, length(curEnGrid)),
  #                                 uncs =  curUncs,
  #                                 opts = list(ens = curEnGrid,
  #                                             order = 3, outsideZero = TRUE))    
}

curSysDtGPobs <- reacHandlerGPobs$createSysDtGpObs()



mapAssignment <- reacHandlerGPobs$getMapAssignment()[,data.table(REAC, EXPID)]

exactGPObsSysDt <- copy(expDt) 
setnames(exactGPObsSysDt, "L1", "EN")
setnames(exactGPObsSysDt, "EXPID", "EXPIDOLD")
exactGPObsSysDt <-exactGPObsSysDt[mapAssignment,on="REAC"]
exactGPObsSysDt[,GPTYPE:="sqrexp"]

reactions <- exactGPObsSysDt[ , .N, by=EXPID][N>1 ,EXPID]

for(curReac in reactions){
  
  meanCs <- exactGPObsSysDt[EXPID==curReac, mean(abs(DATA-DATAREF))-mean(UNC)]
  gpObsHandler$addGP(expid=curReac, sigma=1, len=0.1, nugget=1e-3)
  #gpObsHandler$addGP(expid=curReac, sigma=ifelse(meanCs < 1, 1, meanCs), len=2, nugget=1e-3)
}

gpObsDt <- gpObsHandler$createGPDt()

gpObsDt[,ADJUSTABLE := FALSE]
#gpObsDt[EXPID %in% exactGPObsSysDt[,.N, by=EXPID][N>1,EXPID],  ADJUSTABLE := PARNAME %in% c("sigma","len")]
gpObsDt[EXPID %in% exactGPObsSysDt[,.N, by=EXPID][N>1,EXPID],  ADJUSTABLE := PARNAME %in% c("sigma","len","nugget")]
gpObsDt[EXPID %in% exactGPObsSysDt[,.N, by=EXPID][N==1,EXPID] & PARNAME %in% c("sigma","len"),  PARVAL := 1]
gpObsDt[EXPID %in% exactGPObsSysDt[,.N, by=EXPID][N==1,EXPID] & PARNAME %in% c("nugget"),  PARVAL := 1]
gpObsDt[, IDX := seq_len(.N)]

gpObsHandler$updateSysDt(curSysDtGPobs)
curSysDtGPobs[,IDX:=seq(1,.N)]
curSysDtGPobs[,REFDATA:=rep(0,.N)]
curSysDtGPobs[,ADJUSTABLE:=rep(FALSE,.N)]
curSysDt <- rbind(curSysDt, curSysDtGPobs, fill=TRUE)
curSysDt[, IDX := seq_len(.N)]

###
#Stest <- reacHandlerGPobs$map(expDt, curSysDtGPobs, ret.mat=TRUE)
##
gpDt <- rbind(gpDt, gpObsDt, fill=TRUE)
gpDt[, IDX := seq_len(.N)]

# create global handler and register the individual handlers
sysCompHandler <- createSysCompHandler()
sysCompHandler$addHandler(reacHandlerGPobs)
sysCompHandler$addHandler(normHandler)
sysCompHandler$addHandler(talysHandler)
sysCompHandler$addGPHandler(gpHandler)

# setup optimization specification

optfuns <- createMLOptimFuns()
optfuns$setDts(curExpDt, curSysDt, gpDt,
               sysCompHandler = sysCompHandler)

# define the hyperparameter limits
setkey(gpDt, IDX)
lowerLims <- rep(NA_real_, nrow(gpDt))
upperLims <- rep(NA_real_, nrow(gpDt))

gpDt[grepl("TALYS-",EXPID) & PARNAME == "sigma", LOWERLIMS := 0.1]
gpDt[grepl("TALYS-",EXPID) & PARNAME == "len", LOWERLIMS := 3]
gpDt[grepl("TALYS-",EXPID) & PARNAME == "nugget", LOWERLIMS := 1e-4]

gpDt[grepl("TALYS-",EXPID)  & PARNAME == "sigma", UPPERLIMS := 0.5]
gpDt[grepl("TALYS-",EXPID)  & PARNAME == "len", UPPERLIMS := 5]
gpDt[grepl("TALYS-",EXPID)  & PARNAME == "nugget", UPPERLIMS := 1000]


gpDt[grepl("REACEXP",EXPID) & PARNAME == "sigma" , LOWERLIMS := 1] # Set different upper limit for GPs on the observable
gpDt[grepl("REACEXP",EXPID) & PARNAME == "len" , LOWERLIMS := 0.01]# Set different upper limit for GPs on the observable
gpDt[grepl("REACEXP",EXPID) & PARNAME == "nugget" , LOWERLIMS := 1e-3] # Set different upper limit for GPs on the observable

#gpDt[grepl("REACEXP",EXPID) & PARNAME == "sigma" , UPPERLIMS := 3*PARVAL]
gpDt[grepl("REACEXP",EXPID) & PARNAME == "sigma" , UPPERLIMS := 2000]

gpDt[grepl("REACEXP",EXPID) & PARNAME == "len" , UPPERLIMS := 5] # Set different upper limit for GPs on the observable
gpDt[grepl("REACEXP",EXPID) & PARNAME == "nugget" , UPPERLIMS := 10] # Set different upper limit for GPs on the observable
gpDt[,INITVAL:=PARVAL]

# optimize hyperparameters
# this may take a few minutes
library(optimParallel)
# Setup of multicore optimization using opimparalell
nCores <- detectCores(all.tests = FALSE, logical = TRUE)
cl <- makeCluster(6)
setDefaultCluster(cl=cl)
#optimParallel
optRes <- optimParallel(par = gpDt[ADJUSTABLE==TRUE, INITVAL], 
                        fn = optfuns$logLike, 
                        gr = optfuns$gradLogLike, 
                        method = "L-BFGS-B",
                        lower = gpDt[ADJUSTABLE==TRUE, LOWERLIMS], 
                        upper = gpDt[ADJUSTABLE==TRUE, UPPERLIMS], 
                        control = list(fnscale = -1)
)

#optRes <- optim(par = gpDt[ADJUSTABLE==TRUE, INITVAL], 
#                fn = optfuns$logLike,
#                method = "L-BFGS-B",
#                lower = gpDt[ADJUSTABLE==TRUE, LOWERLIMS], 
#                upper = gpDt[ADJUSTABLE==TRUE, UPPERLIMS], 
#                control = list(fnscale = -1)
#)

newDts <- optfuns$getModifiedDts(optRes$par)

optExpDt <- newDts$expDt
optSysDt <- newDts$sysDt
optGpDt <- newDts$gpDt

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)

