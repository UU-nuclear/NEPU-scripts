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

scriptnr <- 4L
overwrite <- FALSE

library(ggnewscale)
library(reshape2)
library(gplots)
library(RColorBrewer)

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

subents <- read_object(1, "subents")
refParamDt <- read_object(2, "refParamDt")
modList <- read_object(3, "modList")
modDt <- read_object(3, "modDt") # default model prediction mapped to experimental energies
expDt <- read_object(3, "expDt")
sysUncDt <- read_object(3, "sysUncDt")

expDt <- expDt[!is.na(UNC)]

# define objects to be returned
outputObjectNames <- c("fake_expDt", "full_covMat", "fake_subents")
check_output_objects(scriptnr, outputObjectNames)

# try to change the data set that has an absolute systematic uncertainty to a relative one
sysUncDt[ERRTYPE=="sys-abs"]$UNC <- 0.05
sysUncDt[ERRTYPE=="sys-abs"]$ERRTYPE <- "sys-rel"
#sysUncDt[ERRTYPE=="sys-abs"]$UNC <- 10

# instantiate handlers 
reacHandler <- createSysCompReacHandler(c(subents, modList$SUBENT))
normHandler <- createSysCompNormHandler("DATAREF")

# define normalization error of experiments retrieved from EXFOR
# and with some rule based correction/penalization
# (sysUncDt$UNC will be enlarged in this step if necessary using MLO)
normHandler$addSysUnc("EXPID", sysUncDt$EXPID, 
                      0, sysUncDt$UNC, rel = (sysUncDt$ERRTYPE == "sys-rel"))

# define the second derivative prior on the reaction cross sections
reacHandler$addMap("pw", pwMap)

# We define for every reaction channel a prior 
# that incorporates the knowledge that cross sections
# are smooth functions. We achieve this by defining an
# independent Gaussian prior on the second derivatives of the 
# cross section. Because thresholds are different for different 
# reactions, we locate the threshold and cut away the energy  
# grid points below that value.
# Also the prior on the 2nd derivative of the cross section
# in 'curUncs' is modified depending on the reaction channel.
# As a rule of thumb: Reaction channels with larger cross sections
# should permit bigger changes in the 2nd derivative of the cross
# section.
# Meaning of elements in curUncs:
# 1st element: prior uncertainty of value at threshold
# 2nd element: prior uncertainty of the slope of the xs at threshold
# Other elements: prior uncertainty of the 2nd derivative of the xs

tot_length <- 0
EnGridPolyChains <- vector()
for(curReac in unique(expDt$REAC)) {
	# define a spline interpolation of the talys model
	curModDt <- modDt[REAC==curReac]
	spline_func <- splinefun(curModDt$L1,curModDt$DATA)

	threshold_idx <- which(curModDt[,DATA]==0)
	if(length(threshold_idx)==0) {
		#threshold_energy <- curModDt$L1[1] - 0.2
		threshold_idx <- 1
	} else {
		#threshold_energy <- curModDt$L1[max(threshold_idx)]
		threshold_idx <- max(threshold_idx)
	}
	
	#En_min <- min(expDt[REAC==curReac]$L1)
	#threshold_idx <- max(which(curModDt$L1<En_min))

	En_max <- max(expDt[REAC==curReac]$L1)
	max_idx <- which(curModDt$L1>En_max)[1]

	curEnGrid <- curModDt$L1[threshold_idx:max_idx]

	#cat(curReac,"\n")
	#print(threshold_idx)


	#second_derivatives <- spline_func(curEnGrid[1:(length(curEnGrid)-2)],deriv=2)
	y <- spline_func(curEnGrid)
	first_derivatives <- y[2:length(y)] - y[1:(length(y)-1)]
	
	par_vals <- c(y[1],first_derivatives[1])
	if(length(curEnGrid)>2) {
		second_derivatives <- first_derivatives[2:length(first_derivatives)] - first_derivatives[1:(length(first_derivatives)-1)]
		par_vals <- c(y[1],first_derivatives[1],second_derivatives)
	}

	unc_multiplier <- 0.5
	par_uncs <- 0.1*unc_multiplier*pmax(abs(par_vals),1)
	par_uncs[1] <- max(abs(par_vals[1]),0.01)
	par_uncs[2] <- 3*max(abs(first_derivatives[1]),0.01)

	if(length(curEnGrid)>2) {
		#par_uncs[3:length(par_uncs)] <- unc_multiplier*max(abs(second_derivatives))
		par_uncs[3:length(par_uncs)] <- pmax(unc_multiplier*abs(second_derivatives),0.1)
	}

	reacHandler$assignMapToReac("pw", curReac,
                            vals = par_vals,
                            uncs = par_uncs,
                            opts = list(ens = curEnGrid,
                                        order = 3, outsideZero = TRUE))

	tot_length <- tot_length + length(curEnGrid)
	EnGridPolyChains <- c(EnGridPolyChains,curEnGrid)
}

sysCompHandler <- createSysCompHandler()
sysCompHandler$addHandler(reacHandler)
sysCompHandler$addHandler(normHandler)

sysDt <- sysCompHandler$createSysDt()
sysDt[, REFDATA := 0]  
origSysDt <- copy(sysDt)

# test to set a minimum systematic uncertainty on all data set of 1%
# sysDt[ERRTYPE=="sys-rel" & UNC < 0.01]$UNC <- 0.01

# define for MLO optimization that we want to optimize the
# experimental normalization uncertainties and not the uncertainties
# on the second derivatives of the cross section prior
# nor the statistical uncertainties of the experiments

sysDt[, ADJUSTABLE := FALSE]
expDt[, ADJUSTABLE := FALSE]
sysDt[grepl("^EXPID-", EXPID), ADJUSTABLE := TRUE]
# ----------------------------------------------------------------------------------


# try to set REFDATA to something else
# expDt[, REFDATA := 0]

# start the MLO optimization for the individual channels

optfuns <- createMLOptimFuns()
setkey(sysDt, IDX)
sysDt[, ORIGIDX := IDX]

# store modified uncertainties in updSysDt
updSysDt <- copy(sysDt)

# set the seed for random number generator
# to have reproducible results
set.seed(tuneExpUncSeed)

reac_map_assignment <- reacHandler$getMapAssignment()

for (curReac in unique(expDt$REAC)) {

	cat("CURRENT REACTION: ", curReac, "\n")

	curExpIds <- expDt[REAC == curReac, unique(EXPID)]

	# create expDt only containing experiments of
	# current reaction channel
	curExpDt <- expDt[EXPID %in% curExpIds,]
	curExpDt[, IDX := seq_len(.N)]

	cur_EXPID_reac <- reac_map_assignment[REAC==curReac,EXPID]
	# create sysDt only containing systematic errors of
	# experiments in current reaction channel
	curSysDt <- sysDt[EXPID %in% paste0("EXPID-",curExpIds) | 
	                  #grepl("^REACEXP", EXPID),]   
	                  grepl(cur_EXPID_reac, EXPID),]   
	curSysDt[, IDX := seq_len(.N)]

	# set up the current optimization problem
	# upper limit is 50% for relative uncertainty components and
	# 200 mBarns for absolute uncertainty components
	optfuns$setDts(curExpDt, curSysDt, sysCompHandler = sysCompHandler)
	lowerLims <- curSysDt[ADJUSTABLE == TRUE, UNC]
	#lowerLims <- rep(1E-06,length(curSysDt[ADJUSTABLE == TRUE, UNC])) # test to let MLO reduce the systematic uncs
	upperLims <- curSysDt[ADJUSTABLE == TRUE, ifelse(ERRTYPE=="sys-rel", 0.5, 1000)]
	#upperLims <- rep(2E-06,length(curSysDt[ADJUSTABLE == TRUE, UNC])) # test to let MLO reduce the systematic uncs
	initUncs <- curSysDt[ADJUSTABLE == TRUE, UNC]

	# print logLike associated with reference specification
	# following from rule based approach
	cat("logLike before ML optimization: ", optfuns$logLike(initUncs), "\n")
	# do the bare-bone optimization here
	optRes <- optim(initUncs, optfuns$logLike, optfuns$gradLogLike, 
	                method = "L-BFGS-B", lower = lowerLims, upper = upperLims,
	                control = list(fnscale = -1))  # fnscale -1 to maximize instead of minimize 

	# print logLike after correction of systematic uncertainties
	cat("logLike after ML optimization: ", optRes$value, "\n")

	# save the modified systematic uncertainties
	updSysDt[J(curSysDt[ADJUSTABLE==TRUE, ORIGIDX]),
	         UNC := curSysDt[ADJUSTABLE==TRUE, optRes$par]]
}

# Quick before/after comparison
# updSysDt$UNC - sysDt$UNC


# save the needed files for reference
# save_output_objects(scriptnr, outputObjectNames, overwrite)

# plot the result

# extract the systematic uncertainties mapped to the experimental
# data points for plotting
normHandler2 <- createSysCompNormHandler("DATAREF")
normHandler2$addSysUnc("EXPID", "", 0, 0, TRUE)

sysCompHandler2 <- createSysCompHandler()
sysCompHandler2$addHandler(normHandler2)

S <- sysCompHandler2$map(expDt, origSysDt, ret.mat = TRUE)
origX <- sysCompHandler2$cov(origSysDt, ret.mat = TRUE)
updX <- sysCompHandler2$cov(updSysDt, ret.mat = TRUE)
statUnc <- getDt_UNC(expDt)

origUnc <- sqrt(statUnc^2 + diag(S %*% origX %*% t(S))) 
updUnc <- sqrt(statUnc^2 + diag(S %*% updX %*% t(S)))

setkey(expDt, IDX)
expDt[, ORIGUNC := origUnc]
expDt[, UPDUNC := updUnc]

# extract the GP piece-wise linear
# in the notation of the paper Eq. (59) [left hand side]
# P = P_pwl = sysCompHandler$cov(origSysDt[ERRTYPE=="pw"], ret.mat = TRUE)
# S_mod = S_pwl = sysCompHandler$map(expDt, origSysDt[ERRTYPE=="pw"], ret.mat = TRUE)
# U = sysCompHandler$cov(origSysDt, ret.mat = TRUE)
# S = sysCompHandler$map(expDt, ret.mat = TRUE)
#
# this gives the conditional mean of the GP
# (S_mod %*% P)^T %*% (SUS^T + D) %*% sigma_exp =
# t(S_mod %*% P) %*% mult_invCov_x(getDt_DATA(expDt),Diagonal(x=expDt$UNC^2),S,U)
#

# mapping matrix for linear interpolation of the polygon chain
# the UNC column is the same in the original and the updated
# SysDt[ERRTYPE=="pw"]
S_pwl <- sysCompHandler$map(expDt, origSysDt[ERRTYPE=="pw"], ret.mat = TRUE)
# prior on the second derivative of the polygon chain
P_pwl <- sysCompHandler$cov(origSysDt[ERRTYPE=="pw"], ret.mat = TRUE)

# mapping matrix for systematics to experimental data
S <- sysCompHandler$map(expDt, origSysDt, ret.mat = TRUE)
# Systematic uncertainties of experimental data (before MLO) + 
# prior on the second derivative of the polygon chain
U <- sysCompHandler$cov(origSysDt, ret.mat = TRUE)

# mapping matrix for systematics to experimental data after MLO
S_extra <- sysCompHandler$map(expDt, updSysDt, ret.mat = TRUE)
# Systematic uncertainties of experimental data (before MLO) + 
# prior on the second derivative of the polygon chain
U_extra <- sysCompHandler$cov(updSysDt, ret.mat = TRUE)

prior_polygon_chain_parameters <- sysDt[ERRTYPE=="pw"]$DATA
prior_exp_data <- expDt[,DATAREF] # if the prior on the polygon chain parameters are taken from talys
#prior_exp_data <- rep(0,nrow(expDt)) # if the prior on the polygon chain parameters are set to 0
# posterior mean of the polygon chain parameters given the experimental data
# and their covariance
#polygon_chain_parameters <- t(S_pwl %*% P_pwl) %*% mult_invCov_x(getDt_DATA(expDt),Diagonal(x=expDt$UNC^2),S,U)
polygon_chain_parameters <- prior_polygon_chain_parameters +
														t(S_pwl %*% P_pwl) %*% mult_invCov_x(expDt[,DATA]-prior_exp_data,Diagonal(x=expDt$UNC^2),S,U)
updSysDt[ERRTYPE=="pw", DATA := as.vector(polygon_chain_parameters)]


# posterior mean of the polygon chain parameters given the experimental data
# ignoring the systematic errors in the data
#polygon_chain_parameters <- t(S_pwl %*% P_pwl) %*% mult_invCov_x(getDt_DATA(expDt),Diagonal(x=expDt$UNC^2),S_pwl,P_pwl)
polygon_chain_parameters <- prior_polygon_chain_parameters +
														t(S_pwl %*% P_pwl) %*% mult_invCov_x(expDt[,DATA]-prior_exp_data,Diagonal(x=expDt$UNC^2),S_pwl,P_pwl)
updSysDt[ERRTYPE=="pw", V1 := as.vector(polygon_chain_parameters)]

# posterior mean of the polygon chain parameters given the experimental data
# including the MLO added uncertainties
#polygon_chain_parameters <- t(S_pwl %*% P_pwl) %*% mult_invCov_x(getDt_DATA(expDt),Diagonal(x=expDt$UNC^2),S_extra,U_extra)
polygon_chain_parameters <- prior_polygon_chain_parameters +
														t(S_pwl %*% P_pwl) %*% mult_invCov_x(expDt[,DATA]-prior_exp_data,Diagonal(x=expDt$UNC^2),S_extra,U_extra)

updSysDt[ERRTYPE=="pw", V2 := as.vector(polygon_chain_parameters)]

priorCovMat <- Diagonal(x=(sysDt[ERRTYPE=='pw',UNC])^2)
polygon_chain_covMat <- priorCovMat - mult_xt_invCov_x((S_pwl %*% P_pwl),Diagonal(x=expDt$UNC^2),S_extra,U_extra)

# I could replace this mapping below with something like this
# SS <- sysCompHandler$map(curExpDt,curSysDt,ret.mat=TRUE)
# XS_Poly_chain <- SS %*% (curSysDt[,V2])
# curExpDt should here be a data.table with the energies of the knotpoints in the L1 column

covmat_list <- list()
for (curReac in unique(updSysDt[ERRTYPE=="pw"]$EXPID)) {

	curSysDt <- updSysDt[ERRTYPE=="pw" & EXPID==curReac]
	v <- curSysDt$DATA

	idcs <- curSysDt$IDX
	curCovMat <- polygon_chain_covMat[idcs,idcs]
	cur_priorCovMat <- priorCovMat[idcs,idcs]

	T <- matrix(1,nrow=length(v),ncol=length(v))
	T <- lower.tri(T, diag=TRUE) * T

	R <- matrix(1,nrow=length(v),ncol=length(v))
	R <- lower.tri(R, diag=TRUE) * R
	R[,1] <- 0
	R[1,1] <- 1

	Z <- T %*% R

	y <- Z %*% v
	updSysDt[ERRTYPE=="pw" & EXPID==curReac,XS_POLY_CHAIN:=as.vector(y)]

	v <- curSysDt$V1
	y <- Z %*% v
	updSysDt[ERRTYPE=="pw" & EXPID==curReac,XS_POLY_CHAIN_STAT:=as.vector(y)]

	v <- curSysDt$V2
	y <- Z %*% v
	updSysDt[ERRTYPE=="pw" & EXPID==curReac,XS_POLY_CHAIN_MLO:=as.vector(y)]

	covMat_y <- Z %*% curCovMat %*% t(Z)
	updSysDt[ERRTYPE=="pw" & EXPID==curReac,XS_POLY_CHAIN_MLO_UNC:=sqrt(diag(covMat_y))]
	covmat_list[[curReac]] <- covMat_y

	#the prior
	v <- sysDt[ERRTYPE=="pw" & EXPID==curReac]$DATA
	y <- Z %*% v
	updSysDt[ERRTYPE=="pw" & EXPID==curReac,XS_POLY_CHAIN_PRIOR:=as.vector(y)]


	prior_covMat_y <- Z %*% cur_priorCovMat %*% t(Z)
	updSysDt[ERRTYPE=="pw" & EXPID==curReac,XS_POLY_CHAIN_PRIOR_UNC:=sqrt(diag(prior_covMat_y))]
}

updSysDt[ERRTYPE=="pw",EN_POLY_CHAIN:=EnGridPolyChains]

for(row_idx in seq_len(nrow(reac_map_assignment))) {
	expID <- reac_map_assignment[row_idx,EXPID]
	reac <- reac_map_assignment[row_idx,REAC]
	updSysDt[ERRTYPE=="pw" & EXPID==expID, REAC:=reac]
}

# now create a data.table with "experimental" data

# invent some exfor entry numbers
expIDs <- seq(9) + max(as.numeric(expDt[,EXPID]))
reacs <- expDt[,unique(REAC)]

for(i in seq_len(length(reacs))) {
	channel <- reacs[i]

	cur_SysDt <- updSysDt[REAC==channel]

	idx_col <- seq_len(nrow(cur_SysDt))
	id_col <- rep(expIDs[i],times=nrow(cur_SysDt))
	reac_col <- rep(reacs[i],times=nrow(cur_SysDt))
	L1_col <- cur_SysDt[,EN_POLY_CHAIN]
	L2_col <- rep(0,times=nrow(cur_SysDt))
	L3_col <- rep(0,times=nrow(cur_SysDt))
	data_col <- cur_SysDt[,XS_POLY_CHAIN_MLO]

	channel_expDt <- data.table(IDX=idx_col,
															EXPID=id_col,
															REAC=reac_col,
															L1=L1_col,
															L2=L2_col,
															L3=L3_col,
															DATA=data_col)
}

fake_subents <- list()
for(i in seq_len(length(reacs))) {

	channel <- reacs[i]
	cur_SysDt <- updSysDt[REAC==channel]

	subent <- list()
	subent$ID <- expIDs[i]

	subent_data_table <- data.table(
			"EN" = cur_SysDt[,EN_POLY_CHAIN],
			"DATA" = cur_SysDt[,XS_POLY_CHAIN_MLO],
			"ERR-T" = cur_SysDt[,XS_POLY_CHAIN_MLO_UNC]
			)

	subent$DATA$TABLE <- subent_data_table
	subent$DATA$UNIT$EN <- "MEV"
	subent$DATA$UNIT$DATA <- "MB"
	subent$DATA$UNIT$'ERR-T' <- "MB"

	subent$DATA$DESCR <- c("EN","DATA","ERR-T")
	subent$BIB$REACTION <- reacs[i]

	fake_subents[[i]] <- subent
}

fake_expDt <- exforHandler$extractData(fake_subents, ret.values = TRUE)

# create the full covariance matrix for the fake expDt, from the blocks of each
# reaction channel
full_covMat <- bdiag(covmat_list)
fake_expDt[,TOT_UNC:=sqrt(diag(full_covMat))]

# remove points that do not have data
idx_to_keep <- c()
for(channel in reacs) {
	curEnGrid <- fake_expDt[REAC==channel,L1]
	grid_points_to_keep <- c()
	for(grid_i in seq(from=1,to=(length(curEnGrid)-1))) {
		npoints_in_bin <- sum(expDt[REAC==channel,L1] > curEnGrid[grid_i] & expDt[REAC==channel,L1] <= curEnGrid[grid_i+1])
		if(npoints_in_bin) grid_points_to_keep <- c(grid_points_to_keep,grid_i,grid_i+1)
	}
	grid_points_to_keep <- unique(grid_points_to_keep)
	idx_to_keep <- c(idx_to_keep,fake_expDt[REAC==channel,IDX][1] + grid_points_to_keep - 1)
}

fake_expDt <- fake_expDt[idx_to_keep]
fake_expDt[,IDX:=seq_len(.N)]
fake_expDt[,DIDX:=seq_len(.N)]
full_covMat <- full_covMat[idx_to_keep,idx_to_keep]

save_output_objects(scriptnr, outputObjectNames, overwrite)

# it should be fine to use the same needsDt since the energy range of calculations did not change
# this should also mean that the reference Jacobian should be exactly the same (if not there is something wrong)
# I will need to redo all the scripts following the reference Jacobian calculation, since I'm going away from
# having the representeation of the covariance matrix as a diagonal matrix with the random uncertainties and blocks
# with systematic uncs.

# Is hould have two matrices S,U to use mult_invCov_x():
# S is a mapping matrix with dimension: nrow(expDt) x nrow(sysDt)
# U is a diagonal matrix with dimension: nrow(sysDt) x nrow(sysDt)
# If I just set S to be a null-matrix U shouldn't be mapped
# to the experiments at all
#
# Seems to work, the following tests evaluate all(x1 == x2) : TRUE and A2 == A2 : TRUE
# 
# Stest <- Matrix(data=0,nrow=nrow(fake_expDt),ncol=nrow(sysDt))
#
# x1 <- mult_invCov_x(fake_expDt[,DATA],full_covMat,Stest,U)
# x2 <- solve(full_covMat,fake_expDt[,DATA])
#
# A1 <- mult_xt_invCov_x(fake_expDt[,DATA],full_covMat,Stest,U)
# A2 <- t(fake_expDt[,DATA]) %*% x2


# create some plots of the result
ggp1 <- ggplot(expDt) + theme_bw() +
theme(text=element_text(size=2),
		axis.text=element_text(size=2),
	    axis.title=element_text(size=3),
	    plot.title=element_text(size=3),
	    legend.text=element_text(size=2),
	    legend.title=element_text(size=2)) +
xlab("energy [MeV]") + ylab("cross section [mbarn]") +
guides(col = "none") +
geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = "black",
                               linewidth = 0.1, width = 0.3) +
geom_errorbar(aes(x = L1, ymin = DATA - ORIGUNC, ymax = DATA + ORIGUNC, col = EXPID),
                           linewidth = 0.1, width = 0.2) +
geom_point(aes(x = L1, y = DATA, col = EXPID),size=0.2) +

new_scale_colour() +

geom_line(data = updSysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN, col='orig. sys. unc.'), linewidth = 0.1) +
geom_line(data = updSysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN_STAT, col='only stat. unc.'), linewidth = 0.1) +
geom_line(data = updSysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN_MLO, col='extra sys. unc.'), linewidth = 0.1) +
geom_ribbon(data = updSysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, ymin = XS_POLY_CHAIN_MLO - XS_POLY_CHAIN_MLO_UNC,
	ymax = XS_POLY_CHAIN_MLO + XS_POLY_CHAIN_MLO_UNC), fill='green', alpha=0.3) +
#geom_line(data = updSysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN_PRIOR, col='prior'), linewidth = 0.1) +
scale_color_manual(name='Regression Model',
                     breaks=c('orig. sys. unc.', 'only stat. unc.', 'extra sys. unc.','prior'),
                     values=c('extra sys. unc.'='green', 'only stat. unc.'='red', 'orig. sys. unc.'='black','prior'='blue')) +
#ylim(0,NA) +
facet_wrap(~ REAC, ncol=3, scales='free')

dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
filepath <- file.path(plotPath, 'DataReplacement.png')
#ggsave(filepath, ggp1, width = 53, height = 29.8125, units = "cm", dpi = 300)
ggsave(filepath, ggp1, width = 0.2*53, height = 0.2*29.8125, units = "cm", dpi = 300)

plotDir <- file.path(plotPath, 'DataReplacement')
dir.create(plotDir, recursive=TRUE, showWarnings=FALSE)
for(channel in unique(expDt[,REAC])) {
	cur_expDt <- expDt[REAC==channel]
	cur_SysDt <- updSysDt[REAC==channel]
	cur_modDt <- modDt[REAC==channel]

	ggp1 <- ggplot(cur_expDt) + theme_bw() +
	theme(text=element_text(size=2),
			axis.text=element_text(size=2),
		    axis.title=element_text(size=3),
		    plot.title=element_text(size=3),
		    legend.text=element_text(size=2),
		    legend.title=element_text(size=2)) +
	xlab("energy [MeV]") + ylab("cross section [mbarn]") +
	guides(col = "none") +
	geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = "black",
	                               linewidth = 0.1, width = 0.3) +
	geom_errorbar(aes(x = L1, ymin = DATA - ORIGUNC, ymax = DATA + ORIGUNC, col = EXPID),
	                           linewidth = 0.1, width = 0.2) +
	geom_point(aes(x = L1, y = DATA, col = EXPID),size=0.2) +

	new_scale_colour() +
	geom_line(data = cur_SysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN, col='orig. sys. unc.'), linewidth = 0.1) +
	geom_line(data = cur_SysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN_STAT, col='only stat. unc.'), linewidth = 0.1) +
	geom_line(data = cur_SysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN_MLO, col='extra sys. unc.'), linewidth = 0.1) +
	#geom_line(data = cur_modDt, aes(x = L1, y = DATA), col='cyan', linewidth = 0.1) +
	geom_ribbon(data = cur_SysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, ymin = XS_POLY_CHAIN_MLO - XS_POLY_CHAIN_MLO_UNC,
	ymax = XS_POLY_CHAIN_MLO + XS_POLY_CHAIN_MLO_UNC), fill="green", alpha=0.3) +
	geom_line(data = cur_SysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN_PRIOR, col='prior'), linewidth = 0.1) +
	#geom_ribbon(data = cur_SysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, ymin = XS_POLY_CHAIN_PRIOR - XS_POLY_CHAIN_PRIOR_UNC,
	#ymax = XS_POLY_CHAIN_PRIOR + XS_POLY_CHAIN_PRIOR_UNC), fill="blue", alpha=0.3) +
	
	#geom_line(data = cur_SysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN_PRIOR), linewidth = 0.1, col='cyan') +
	scale_color_manual(name='Regression Model',
	                     breaks=c('orig. sys. unc.', 'only stat. unc.', 'extra sys. unc.','prior'),
	                     values=c('extra sys. unc.'='green', 'only stat. unc.'='red', 'orig. sys. unc.'='black','prior'='blue'))



	filename <- file.path(plotDir, paste0(channel,'.png'))
	ggsave(filename, ggp1, width = 0.2*53, height = 0.2*29.8125, units = "cm", dpi = 300)

	expID <- reac_map_assignment[REAC==channel,EXPID]
	curCovMat <- covmat_list[[expID]]

	# save plots of the correlation matrices also
	ncolors <- 256
	color_palette<-colorRampPalette(c("red","white","blue"))(ncolors)
	png(file=file.path(plotDir, paste0(channel,'_corrMat.png')),
	width=1200, height=700)
	heatmap.2(cov2cor(as.matrix(curCovMat)),
	  Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",
	  margins=c(8,8),trace = "none",dendrogram="none",density.info = "none",
	  breaks=seq(-1.,1.,length.out = ncolors+1))
	dev.off()
}

# save a plot of the final full correlation matrices also
ncolors <- 256
color_palette<-colorRampPalette(c("red","white","blue"))(ncolors)
png(file=file.path(plotDir, paste0('Full_corrMat.png')),
width=1200, height=700)
heatmap.2(cov2cor(as.matrix(full_covMat)),
  Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",
  margins=c(8,8),trace = "none",dendrogram="none",density.info = "none",
  breaks=seq(-1.,1.,length.out = ncolors+1))
dev.off()

ggp1 <- ggplot(expDt) + theme_bw() +
theme(text=element_text(size=2),
		axis.text=element_text(size=2),
	    axis.title=element_text(size=3),
	    plot.title=element_text(size=3),
	    legend.text=element_text(size=2),
	    legend.title=element_text(size=2)) +
xlab("energy [MeV]") + ylab("cross section [mbarn]") +
guides(col = "none") +
geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = "black",
                               linewidth = 0.1, width = 0.3) +
geom_errorbar(aes(x = L1, ymin = DATA - ORIGUNC, ymax = DATA + ORIGUNC, col = EXPID),
                           linewidth = 0.1, width = 0.2) +
geom_point(aes(x = L1, y = DATA, col = EXPID),size=0.2) +
geom_ribbon(data = fake_expDt,aes(x = L1, ymin = DATA - TOT_UNC,
	ymax = DATA + TOT_UNC), fill='green', alpha=0.3) +
geom_errorbar(data = fake_expDt,aes(x = L1, ymin = DATA - TOT_UNC,
	ymax = DATA + TOT_UNC), col='black',
                           linewidth = 0.1, width = 0.2) +
facet_wrap(~ REAC, ncol=3, scales='free')