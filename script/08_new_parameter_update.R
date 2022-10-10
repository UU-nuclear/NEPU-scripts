#
# DESCRIPTION OF STEP
#

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

scriptnr <- 8L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

refParamDt <- read_object(2, "refParamDt")
optGpDt <- read_object(6, "optGpDt")
optParamDt <- read_object(7, "optParamDt")
yexp <- read_object(7, "yexp")
D <- read_object(7, "D")
S0 <- read_object(7, "S0")
X <- read_object(7, "X")
optRes <- read_object(7, "optRes")
optSysDt_allpars <- read_object(7, "optSysDt_allpars")
P0 <- read_object(7, "P0")

##################################################
#       START OF SCRIPT
##################################################
print("-----------------------------------------------------")
print("----------------------script 08----------------------")
print("-----------------------------------------------------")

# define objects to be returned
outputObjectNames <- c("finalPars_new", "finalParamDt_new")
check_output_objects(scriptnr, outputObjectNames)

# construct the full prior covariance matrix
# including sensitive and insensitive parameters
# sensitive parameters were considered during LM optimization
# whereas insensitive parameters not
gpHandler <- createSysCompGPHandler()
sysCompHandler <- createSysCompHandler()
sysCompHandler$addGPHandler(gpHandler)
P0_all <- sysCompHandler$cov(optSysDt_allpars, optGpDt, ret.mat = TRUE)

# update insensitive parameters
# not considered for optimization

p0 <- unlist(refParamDt[ADJUSTABLE==TRUE, PARVAL])

# p_upd is the vector of updated paramters
d1 <- yexp - optRes$fn
delta_p <- mult_invCov_x(d1, D, S0, X)
#invP0_all <- solve(P0_all)
Jacobian_opt_pars <- optRes$jac
p_upd <- delta_p + p0
finalPars <- p_upd
finalPars[optpars_indices] <- optRes$par


# IMPORTANT REMARK: 
# If Levenberg-Marquardt did not converge sufficiently 
# to a (local) posterior maximum, some values in 
# finalParCovmat may be negative. Such a result does
# not correspond to a proper covariance matrix.

# save the result also in a more human-readable way
# as a table also containing the parameter names
finalParamDt_new <- copy(optParamDt)
finalParamDt_new[PARNAME %in% optSysDt_allpars$PARNAME, ADJUSTABLE:=TRUE]
setkey(finalParamDt_new, IDX)
finalParamDt_new[ADJUSTABLE == TRUE, POSTVAL := paramTrafo$fun(finalPars)]
# IMPORTANT NOTE: POSTUNC is still with respect to transformed parameters
finalParamDt_new[ADJUSTABLE == TRUE, POSTUNC := sqrt(diag(finalParCovmat))]

# save the needed files for reference
# save_output_objects(scriptnr, outputObjectNames, overwrite)
