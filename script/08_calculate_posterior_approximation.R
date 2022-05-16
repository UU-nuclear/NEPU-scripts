#
# DESCRIPTION OF STEP
#
# Approximate the posterior distribution
# by constructing a second order Taylor
# polynomial at the parameter set 
# found by the Levenberg-Marquardt algorithm.
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

scriptnr <- 8L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

refParamDt <- read_object(2, "refParamDt")
optSysDt <- read_object(6, "optSysDt")
optGpDt <- read_object(6, "optGpDt")
extNeedsDt <- read_object(2, "extNeedsDt")
optParamDt <- read_object(7, "optParamDt")
Sexp <- read_object(7, "Sexp")
mask <- read_object(7, "mask")
refPar <- read_object(7, "refPar")
P0 <- read_object(7, "P0")
yexp <- read_object(7, "yexp")
D <- read_object(7, "D")
S0 <- read_object(7, "S0")
X <- read_object(7, "X")
optRes <- read_object(7, "optRes")
optSysDt_allpars <- read_object(7, "optSysDt_allpars")
optSysDt_optpars <- read_object(7, "optSysDt_optpars")

##################################################
#       START OF SCRIPT
##################################################

# define objects to be returned
outputObjectNames <- c("variationMat", "variationMatRes", 
                       "finalPars", "finalParCovmat")
check_output_objects(scriptnr, outputObjectNames)

# see step 07_tune_talyspars.R for more explanation
# about setting up the talys handler
talysHnds <- createTalysHandlers()
talys <- talysHnds$talysOptHnd
talys$setPars(optParamDt)
talys$setParTrafo(paramTrafo$fun, paramTrafo$jac)
talys$setNeeds(extNeedsDt)
talys$setSexp(Sexp)
talys$setMask(mask)
talys$setEps(0.01)

# construct an object that contains functions to
# evaluate the log posterior density (up to a constant)
# and its gradient. It makes use of the talys handler
# and takes in addition the parameters defining the
# posterior distribution.
logPost <- setupLogPost(refPar, P0, yexp, D, S0, X, talys)

# construct the full prior covariance matrix
# including sensitive and insensitive parameters
# sensitive parameters were considered during LM optimization
# whereas insensitive parameters not
gpHandler <- createSysCompGPHandler()
sysCompHandler <- createSysCompHandler()
sysCompHandler$addGPHandler(gpHandler)
P0_all <- sysCompHandler$cov(optSysDt_allpars, optGpDt, ret.mat = TRUE)

# get the indices of the adjustable parameters
optpars_indices <- optSysDt_optpars[, sort(IDX)]

# consistency check
stopifnot(all(P0_all[optpars_indices, optpars_indices] == P0))

# compute the GLS Hessian
# Smod maps from model parameters to optimzed parameters only
Smod <- optRes$jac
tS_invCexp_S <- mult_xt_invCov_x(Smod, D, S0, X)
invP0_all <- solve(P0_all)
H_gls <- (-solve(P0_all))
H_gls[optpars_indices, optpars_indices] <- H_gls[optpars_indices, optpars_indices] - (tS_invCexp_S)

# define the variations needed to calculate the diagonal
# elements of the Hessian matrix
workEps <- 0.01
optPars <- optRes$par

idx1 <- seq_along(optPars) 
idx2 <- 1+seq_along(optPars)
idx3 <- 1+length(optPars)+seq_along(optPars)

variationMat <- matrix(optPars, nrow = length(optPars), ncol = 1+length(optPars)*2)
variationMat[cbind(idx1,idx2)] <- variationMat[cbind(idx1,idx2)] - workEps
variationMat[cbind(idx1,idx3)] <- variationMat[cbind(idx1,idx3)] + workEps

# perform the calculations using the varied parameter sets
variationMatRes <- talys$fun(variationMat)
logPostVals <- logPost$fun(variationMat, variationMatRes)

# calculate the diagonal of the Hessian matrix
# by evaluating the 2nd derivative using finite differences
fp <- (logPostVals[idx3] - logPostVals[idx2]) / (2*workEps)
fpp <- (logPostVals[idx3] - 2 * logPostVals[1] + logPostVals[idx2]) / workEps^2

# compute the approximation to the Hessian
H <- H_gls
diag(H)[optpars_indices] <- pmin(fpp, diag(H_gls)[optpars_indices])

# The second derivative (Hessian matrix) of the posterior distribution
# is approximately the negative inverse covariance matrix.
finalParCovmat <- (-1) * solve(H)

# if LM algorithm did not sufficiently converge
# in step 07, the covariance matrix would not be well defined
stopifnot(isSymmetric(finalParCovmat, tol=1e-9))
finalParCovmat <- (finalParCovmat + t(finalParCovmat)) / 2
stopifnot(all(eigen(finalParCovmat)$values >= 0))

# update insensitive parameters
# not considered for optimization
setkey(refParamDt, IDX)
p0 <- unlist(refParamDt[ADJUSTABLE==TRUE, PARVAL])
pref <- p0
pref[optpars_indices] <- optRes$par

d0 <- pref - p0
d1 <- yexp - optRes$fn

invCexp_d1 <- mult_invCov_x(d1, D, S0, X)
G <- (-invP0_all %*% d0)
G[optpars_indices] <- G[optpars_indices] + t(Smod) %*% invCexp_d1

p1 <- as.matrix(pref + finalParCovmat %*% G)

finalPars <- p1
finalPars[optpars_indices] <- optRes$par

# IMPORTANT REMARK: 
# If Levenberg-Marquardt did not converge sufficiently 
# to a (local) posterior maximum, some values in 
# finalParCovmat may be negative. Such a result does
# not correspond to a proper covariance matrix.

# save the result also in a more human-readable way
# as a table also containing the parameter names
finalParamDt <- copy(optParamDt)
finalParamDt[PARNAME %in% optSysDt_allpars$PARNAME, ADJUSTABLE:=TRUE]
setkey(finalParamDt, IDX)
finalParamDt[ADJUSTABLE == TRUE, POSTVAL := paramTrafo$fun(finalPars)]
# IMPORTANT NOTE: POSTUNC is still with respect to transformed parameters
finalParamDt[ADJUSTABLE == TRUE, POSTUNC := sqrt(diag(finalParCovmat))]

# sanity check:
# make sure that the results obtained by the two different
# talys handlers talysOptHnd and talysHnd coincide.
# This is indication that we have indeed used the correct parameters
# in the optimization routine. In principle, it is still possible
# that the matrix Sexp is not correctly specified. This can be
# checked by looking at the plots comparing predictions and
# experimental data.
# finalParamDt[!is.na(POSTVAL), PARVAL := as.list(POSTVAL)]
# myInp <- convertToInput(finalParamDt, extNeedsDt)
# res1 <- talysHnds$talysOptHnd$fun(finalPars)
# runObj <- talysHnds$talysHnd$run(myInp$inputs, myInp$outspecs)
# talysHnds$talysHnd$isRunning(runObj)
# rawRes2 <- talysHnds$talysHnd$result(runObj)
# res2 <- as.vector(Sexp %*% rawRes2[[1]]$result$V1)
# stopifnot(all(res1 == res2))
# sanity check end

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)
