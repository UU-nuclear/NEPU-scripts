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

if(length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
} else if (length(args) > 1) {
  stop("Script only accepts one argument.", call.=FALSE)
}

#################################################
#       SCRIPT PARAMETERS
##################################################

scriptnr <- 11L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

refParamDt <- read_object(2, "refParamDt")
#optSysDt <- read_object(6, "optSysDt")
optSysDt <- read_object(4, "updSysDt")
#optGpDt <- read_object(6, "optGpDt")
extNeedsDt <- read_object(2, "extNeedsDt")
#optParamDt <- read_object(10, "optParamDt")
optParamDt <- read_object(5, "optParamDt")
Sexp <- read_object(5, "Sexp")
mask <- read_object(5, "mask")
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
outputObjectNames <- c("finalPars", "finalParCovmat","finalParamDt")
check_output_objects(scriptnr, outputObjectNames)

# construct the full prior covariance matrix
# including sensitive and insensitive parameters
# sensitive parameters were considered during LM optimization
# whereas insensitive parameters not
gpHandler <- createSysCompGPHandler()
sysCompHandler <- createSysCompHandler()
sysCompHandler$addGPHandler(gpHandler)
#P0_all <- sysCompHandler$cov(optSysDt_allpars, optGpDt, ret.mat = TRUE)
P0_all <- sysCompHandler$cov(optSysDt_allpars, ret.mat = TRUE)

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

# The second derivative (Hessian matrix) of the posterior distribution
# is approximately the negative inverse covariance matrix.
finalParCovmat <- (-1) * solve(H_gls)

# if LM algorithm did not sufficiently converge
# in step 07, the covariance matrix would not be well defined
stopifnot(isSymmetric(finalParCovmat, tol=1e-8))
finalParCovmat <- (finalParCovmat + t(finalParCovmat)) / 2
stopifnot(all(eigen(finalParCovmat)$values >= 0))

setkey(refParamDt, IDX)
p0 <- unlist(refParamDt[ADJUSTABLE==TRUE, PARVAL])

# prior expectation on optimizied parameters
p0_opt <- p0[optpars_indices]
# optimized parameter vector
p_opt <- optRes$par
# prior expectation on insensitive parameters
p0_insensitive <- p0[-optpars_indices]
# covariance matrix of the optimized parameters
Sigma_22 <- finalParCovmat[optpars_indices,optpars_indices]
# (prior) covariance between insensitive and optimized parameters
Sigma_12 <- finalParCovmat[-optpars_indices,optpars_indices]
# updated insensitive paramters according to the conditional MVN
p_insensitive <- p0_insensitive + Sigma_12 %*% solve(Sigma_22,p_opt - p0_opt)

finalPars <- unlist(refParamDt[ADJUSTABLE==TRUE,PARVAL])
finalPars[optpars_indices] <- p_opt
finalPars[-optpars_indices] <- p_insensitive

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
# We need another parameter transformation here since the full parameter
# vector is different from the optimized parameter vector
adjParNamesFull <- finalParamDt[ADJUSTABLE==TRUE,PARNAME]
for( i in 1:nrow(parRanges) ) {
  finalParamDt[grepl(parRanges[i]$keyword,PARNAME),PARMIN:=parRanges[i]$min]
  finalParamDt[grepl(parRanges[i]$keyword,PARNAME),PARMAX:=parRanges[i]$max]
}

# set the parameter transformation to be centered at the prior mean/mode
# the ranges of the parameters are the ones specified in the TALYS manual
paramTrafo <- parameterTransform(
                  x0 = unlist(finalParamDt[ADJUSTABLE==TRUE,PARVAL]),
                  delta = finalParamDt[ADJUSTABLE==TRUE,unlist(PARVAL) - PARMIN])


finalParamDt[ADJUSTABLE == TRUE, POSTVAL := paramTrafo$fun(finalPars)]
# IMPORTANT NOTE: POSTUNC is still with respect to transformed parameters
#finalParamDt[ADJUSTABLE == TRUE, POSTUNC := sqrt(diag(finalParCovmat))]
par_unc_int <- sqrt(diag(finalParCovmat))
par_min_int <- finalPars-par_unc_int
par_max_int <- finalPars+par_unc_int
par_min_ext <- paramTrafo$fun(par_min_int)
par_max_ext <- paramTrafo$fun(par_max_int)

finalParamDt[ADJUSTABLE==TRUE,POSTUNC_Low := POSTVAL - par_min_ext]
finalParamDt[ADJUSTABLE==TRUE,POSTUNC_UP := par_max_ext - POSTVAL]

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)
