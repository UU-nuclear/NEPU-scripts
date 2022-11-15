#
# DESCRIPTION OF STEP
#
# tune the hyperparameters of the Gaussian
# processes attached to the energy-dependent
# TALYS parameters using maximunm likelihood
# optimization.
#
# This script is the same as 06_tune_endep_hyperpars_fixed_sigma.R
# with the difference that the upper limit for sigma is set by the prior
# uncertainty of the corresponding energy independent parameter

#################################################
#       SCRIPT Setup
##################################################

args = commandArgs(trailingOnly=TRUE)

if (length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}

#################################################
#       SCRIPT PARAMETERS
##################################################

scriptnr <- 6L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

subents <- read_object(1, "subents")
refParamDt <- read_object(2, "refParamDt")
extNeedsDt <- read_object(2, "extNeedsDt")
modList <- read_object(3, "modList")
expDt <- read_object(3, "expDt")
updSysDt <- read_object(4, "updSysDt")
fullSensDt <- read_object(5, "fullSensDt") 

source("misc-scripts/MLOwithPrior.R")
library(parallel)
##################################################
#       START OF SCRIPT
##################################################



print("-----------------------------------------------------")
print("----------------------script 06----------------------")
print("-----------------------------------------------------")

# define objects to be returned
outputObjectNames <- c("optExpDt", "optSysDt", "optGpDt","sysCompHandler")
check_output_objects(scriptnr, outputObjectNames)

############# First: check which parameters are sensitive to data ###########
# convert the sparse matrix given as data.table 
# into a spase matrix type as defined in package Matrix
Spar <- with(fullSensDt,
             sparseMatrix(i = IDX1, j = IDX2, x = X,
                          dims = c(nrow(extNeedsDt), nrow(refParamDt))))
Sexp <- exforHandler$getJac(expDt, extNeedsDt, subents)
Sglob <- Sexp %*% Spar 

# safeguard
stopifnot(all(dim(Sglob) == c(nrow(expDt), nrow(refParamDt))))

# convert the sparse matrix Sglob into a datatable
SglobDt <- as.data.table(summary(Sglob))
setnames(SglobDt, c("IDX1", "IDX2", "X"))

# we (linearly) propagate all parameter values equal one
# to the model predictions
imp1 <- as.vector(Spar %*% rep(1, nrow(refParamDt)))
# we propagate hypothetical experimental values equal one
# to the model prediction
imp2 <- as.vector(t(Sexp) %*% rep(1, nrow(Sexp)))
# then we select observables on the model grid that
# are affected by both the backpropagation from the
# experiment and the forward propagation of model parameters
impIdx <- which(imp1 * imp2 != 0)

optSparDt <- copy(fullSensDt)
setkey(optSparDt, IDX1)
optSparDt <- optSparDt[J(impIdx)]

paramImpactDt <- SglobDt[, list(IMP = max(abs(X))), by = "IDX2"]
paramImpactDt <- paramImpactDt[order(IMP, decreasing = TRUE)]
selParIdcs <- paramImpactDt[IMP >= 1, IDX2]

setkey(optSparDt, IDX2)
mask <- optSparDt[J(selParIdcs), list(DSTIDX = IDX1, SRCIDX = IDX2)]
setkey(mask, SRCIDX, DSTIDX)
adjParIdcs <- unique(mask$SRCIDX)

# make a copy of the reference parameter datatable 
# and define what and what not we want to optimize
optParamDt <- copy(refParamDt)
setkey(optParamDt, IDX)
optParamDt[, ADJUSTABLE := FALSE]
optParamDt[J(adjParIdcs), ADJUSTABLE := TRUE]

# safeguard
stopifnot(sum(optParamDt$ADJUSTABLE) == length(adjParIdcs))

# find all energy dependent parameters which have at least one point that is adjustbale
adjustable_endep_par_names <- optParamDt[ADJUSTABLE==TRUE]$PARNAME[grepl("\\(.+\\)",optParamDt[ADJUSTABLE==TRUE]$PARNAME)]
adjustable_endep_par_names <- unique(str_remove(adjustable_endep_par_names,"\\(.+\\)"))
optParamDt$tmp = str_remove(optParamDt$PARNAME,"\\(.+\\)")
optParamDt[tmp %in% adjustable_endep_par_names]$ADJUSTABLE=TRUE
optParamDt[,tmp:=NULL] # remove the temporary column from the data table


#############################################################################

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
    gpHandler$addGP(curParname, curSysDt[EXPID==curParname,UNC][1], 100, 1e-4)
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
gpDt[, ADJUSTABLE := PARNAME %in% c("len","sigma")]
gpDt[, IDX := seq_len(.N)]

# set hyper-parameters for energy dependent parameters that are not adjustable
# to not be part of the hyper-parameter optiization
gpDt[!(str_remove(EXPID,"TALYS-") %in% adjustable_endep_par_names), ADJUSTABLE:=FALSE]

# initialize handler for systematic experimental uncertainties
# TODO: It should not be needed to define the systematic error
#       components during the setup. The category (column name)
#       should be enough.
normHandler <- createSysCompNormHandler("DATAREF")
tmpDt <- updSysDt[grepl("EXPID-", EXPID)]
expIds <- sub("EXPID-", "", tmpDt$EXPID)
normHandler$addSysUnc("EXPID", expIds, 0, tmpDt$UNC,
                      rel = (tmpDt$ERRTYPE == "sys-rel"))

# create global handler and register the individual handlers
sysCompHandler <- createSysCompHandler()
sysCompHandler$addHandler(normHandler)
sysCompHandler$addHandler(talysHandler)
sysCompHandler$addGPHandler(gpHandler)

# setup optimization specification

optfuns <- createMLOptimFuns()
optfuns$setDts(curExpDt, curSysDt, gpDt,
               sysCompHandler = sysCompHandler)

# create a tempory data table with the uncertainty for each parameter
tmp <- curSysDt[ERRTYPE=="talyspar_endep"]
tmp[,DIDX:=NULL]
tmp[,IDX:=NULL]
tmp[,EN:=NULL]
tmp[,ORIGIDX:=NULL]
tmp <- unique(tmp)
stopifnot(gpDt[PARNAME=="sigma"]$EXPID == tmp$EXPID)
gpDt[PARNAME=="sigma",PARUNC:=tmp$UNC]
# They should have the same order, but should make some reordering of tmp



# define the hyperparameter limits
setkey(gpDt, IDX)
lowerLims <- rep(NA_real_, nrow(gpDt))
upperLims <- rep(NA_real_, nrow(gpDt))
lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- 0.1*gpDt[ADJUSTABLE==TRUE & PARNAME=="sigma",PARUNC]
lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- 3
upperLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- gpDt[ADJUSTABLE==TRUE & PARNAME=="sigma",PARUNC]
upperLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- 200

# prior on the hyper-parameter
priorExpectation <- rep(NA_real_, nrow(gpDt))
priorExpectation[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- gpDt[ADJUSTABLE==TRUE & PARNAME=="sigma",PARUNC]
priorExpectation[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- 50

priorUnc <- rep(NA_real_, nrow(gpDt))
priorUnc[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- 0.5*gpDt[ADJUSTABLE==TRUE & PARNAME=="sigma",PARUNC]
priorUnc[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- 50

#########################

lowerLims <- lowerLims[gpDt$ADJUSTABLE]
upperLims <- upperLims[gpDt$ADJUSTABLE]
priorExpectation <- priorExpectation[gpDt$ADJUSTABLE]
priorUnc <- priorUnc[gpDt$ADJUSTABLE]
initPars <- priorExpectation

# we assume no prior correlations
priorCovMat <- Diagonal(length(priorUnc),x=priorUnc^2)

# set the prior
optfuns$setPrior(priorExpectation,priorCovMat)


library(optimParallel)
# Setup of multicore optimization using optimparalell
nCores <- detectCores(all.tests = FALSE, logical = TRUE)
cl <- makeCluster(min(nCores,32))
setDefaultCluster(cl=cl)
clusterEvalQ(cl, c(library(data.table)))
clusterExport(cl, c("optfuns","logLike","gradLogLike"), 
              envir=environment())
#optimParallel
cat("Started optimization at", as.character(Sys.time()), "\n") 
optRes <- optimParallel(par = initPars, 
                        fn = optfuns$logPost, 
                        gr = optfuns$gradLogpost, 
                        method = "L-BFGS-B",
                        lower = lowerLims, 
                        upper = upperLims, 
                        control = list(fnscale = -1,maxit=300)
)
cat("Finished optimization at", as.character(Sys.time()), "\n")

newDts <- optfuns$getModifiedDts(optRes$par)

optExpDt <- newDts$expDt
optSysDt <- newDts$sysDt
optGpDt <- newDts$gpDt
optGpDt <- optGpDt[order(IDX)]
optGpDt[ADJUSTABLE==TRUE,INITVAL:=initPars]
 
# # sample starting values
# n_sigmas <- length(lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"])
# n_lengths <- length(lowerLims[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"])
# 
# n_adjustable_hyper_pars <- n_sigmas + n_lengths
# hyper_par_vector <- rep(0, n_adjustable_hyper_pars)
# 
# n_samples <- 32
# sigma_samples <- runif(n_sigmas*n_samples,min_sigma,max_sigma)
# length_samples <- runif(n_lengths*n_samples,min_length,max_length)
# 
# hyper_par_samples <- rep(NA_real_, nrow(gpDt)*n_samples)
# hyper_par_samples[gpDt$ADJUSTABLE & gpDt$PARNAME == "sigma"] <- sigma_samples
# hyper_par_samples[gpDt$ADJUSTABLE & gpDt$PARNAME == "len"] <- length_samples
# hyper_par_samples <- hyper_par_samples[!is.na(hyper_par_samples)] # removing the ones that are not adjustable
# hyper_par_samples <- matrix(hyper_par_samples,ncol=n_samples)
# 
# cat("number of hyper-parameters to optimize: ",nrow(gpDt[ADJUSTABLE==TRUE]),"\n")
# 
# nCores <- detectCores(all.tests = FALSE, logical = TRUE)
# cat("number of cores = ",getOption("mc.cores", nCores),"\n")
# cat("Started optimization at", as.character(Sys.time()), "\n") 
# results <- mclapply(X=as.data.frame(hyper_par_samples),FUN=optim, fn=optfuns$logPost, gr=optfuns$gradLogPost,
      # method = "L-BFGS-B", lower = lowerLims, upper = upperLims, control = list(fnscale = -1, maxit=300),
      # mc.cores = getOption("mc.cores", nCores))
# cat("Finished optimization at", as.character(Sys.time()), "\n")
# 
# results <- simplify2array(results)
# best_result <- results[,unlist(results["value",]) == max(unlist(results["value",]))]
# 
# newDts <- optfuns$getModifiedDts(best_result$par)
# #
# optExpDt <- newDts$expDt
# optSysDt <- newDts$sysDt
# optGpDt <- newDts$gpDt

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)

print("---- optimised endep hyperparameters -----")
print(optGpDt[ADJUSTABLE==TRUE],topn=nrow(optGpDt[ADJUSTABLE==TRUE]))
