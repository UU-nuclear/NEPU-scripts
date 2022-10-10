
# data.table containing the prior mean of all parameters
refParamDt <- read_object(2, "refParamDt")

# data.table containing the hyper-paramteres of the energy
# dependent talys parameters
optGpDt <- read_object(6, "optGpDt")

# data.table containing the optimized parameters
optParamDt <- read_object(7, "optParamDt")

# vector of experimental xs
yexp <- read_object(7, "yexp")

# stat. uncertainties
D <- read_object(7, "D")

# mapping of systematic components to experimental data
S0 <- read_object(7, "S0")

# Systematic unc. components
X <- read_object(7, "X")

# data.table containing information about the indexes of optimized parameters
optSysDt_optpars <- read_object(7, "optSysDt_optpars")
optSysDt_allpars <- read_object(7, "optSysDt_allpars")

# results of the LM algorithm optimization
optRes <- read_object(7, "optRes")
#optSysDt_allpars <- read_object(7, "optSysDt_allpars")
Jacobian_opt_pars <- optRes$jac

# construct the full prior covariance matrix
# including sensitive and insensitive parameters
# sensitive parameters were considered during LM optimization
# whereas insensitive parameters not
gpHandler <- createSysCompGPHandler()
sysCompHandler <- createSysCompHandler()
sysCompHandler$addGPHandler(gpHandler)
P0_all <- sysCompHandler$cov(optSysDt_allpars, optGpDt, ret.mat = TRUE)

# p0_all is the prior mean of all paramters, adjustable or not
# ADJUSTABLE==TRUE just deselects the first two rows
# of the data table that includes the energy grid and
# wether or not to produce an endf file
p0_all <- unlist(refParamDt[ADJUSTABLE==TRUE, PARVAL])

# construct the full Jacobian
Jacobian_all_pars <- Matrix(0,nrow=nrow(Jacobian_opt_pars),ncol=length(p0_all))
optpars_indices <- optSysDt_optpars[, sort(IDX)]

# only the entries related to optimized parameters are non-zero
Jacobian_all_pars[,optpars_indices] <- Jacobian_opt_pars

# finally we update the non-sensitive parameters
InvCov_x_residual <- mult_invCov_x(yexp - optRes$fn, D, S0, X)

delta_p <- P0_all %*% t(Jacobian_all_pars) %*% InvCov_x_residual

p_updated <- P0_all %*% t(Jacobian_all_pars) %*% InvCov_x_residual + p0_all

finalPars <- p_updated
finalPars[optpars_indices] <- optRes$par

endep_par_indices <- optSysDt_allpars[ERRTYPE=="talyspar_endep", sort(IDX)]

##############################################

vv <- t(Jacobian_opt_pars) %*% InvCov_x_residual

