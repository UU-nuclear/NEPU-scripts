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

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

subents <- read_object(1, "subents")
modList <- read_object(3, "modList")
modDt <- read_object(3, "modDt") # default model prediction mapped to experimental energies
expDt <- read_object(3, "expDt")
sysUncDt <- read_object(3, "sysUncDt")

expDt <- expDt[!is.na(UNC)]

# define objects to be returned
outputObjectNames <- c("origSysDt", "updSysDt")
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
		threshold_energy <- curModDt$L1[1]
	} else {
		threshold_energy <- curModDt$L1[max(threshold_idx)]
	}

	En_max <- max(expDt[REAC==curReac]$L1) + 0.1
	#curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), En_max, by = 0.1)
	#curEnGrid <- seq(curModDt$L1[1]-0.2, En_max, by = 0.1)
	curEnGrid <- seq(threshold_energy, En_max, by = 0.1)
	#curEnGrid <- seq(0., En_max, by = 0.1)

	# The model tends to fail when using the threshold from talys, i.e getThresEn():
	# at least when there is experimental data below or close to the threshold
	# probably because the accuracy of this threshold depends on the talys evaluation
	# energy grid, which is not very accurate.

	#cat("-------- ",curReac," --------\n")
	#cat("threshold = ", getThresEn(curReac, modDt, defaultThresEn),"\n")
	#cat("min. exp En = ", min(expDt[REAC==curReac,L1]),"\n")

	#second_derivatives <- spline_func(curEnGrid[1:(length(curEnGrid)-2)],deriv=2)
	y <- spline_func(curEnGrid)
	first_derivatives <- y[2:length(y)] - y[1:(length(y)-1)]
	second_derivatives <- first_derivatives[2:length(first_derivatives)] - first_derivatives[1:(length(first_derivatives)-1)]

	#par_vals <- c(spline_func(curEnGrid[1]),spline_func(curEnGrid[1],deriv=1),second_derivatives)
	par_vals <- c(y[1],first_derivatives[1],second_derivatives)
	par_uncs <- par_vals
	par_uncs[1] <- max(10*abs(y[1]),1)
	par_uncs[2] <- max(abs(first_derivatives),1)
	par_uncs[3:length(par_uncs)] <- max(abs(second_derivatives),1)

	#par_vals[2:length(par_uncs)] <- 0

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

for (curReac in unique(expDt$REAC)) {
	cat("CURRENT REACTION: ", curReac, "\n")

	curExpIds <- expDt[REAC == curReac, unique(EXPID)]

	# create expDt only containing experiments of
	# current reaction channel
	curExpDt <- expDt[EXPID %in% curExpIds,]
	curExpDt[, IDX := seq_len(.N)]

	# create sysDt only containing systematic errors of
	# experiments in current reaction channel
	curSysDt <- sysDt[EXPID %in% paste0("EXPID-",curExpIds) | 
	                  grepl("^REACEXP", EXPID),]   
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
save_output_objects(scriptnr, outputObjectNames, overwrite)

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
# posterior mean of the polygon chain parameters given the experimental data
# and their covariance
polygon_chain_parameters <- prior_polygon_chain_parameters +
														t(S_pwl %*% P_pwl) %*% mult_invCov_x(expDt[,DATA]-prior_exp_data,Diagonal(x=expDt$UNC^2),S,U)
updSysDt[ERRTYPE=="pw", DATA := as.vector(polygon_chain_parameters)]

# posterior mean of the polygon chain parameters given the experimental data
# ignoring the systematic errors in the data
polygon_chain_parameters <- prior_polygon_chain_parameters +
														t(S_pwl %*% P_pwl) %*% mult_invCov_x(expDt[,DATA]-prior_exp_data,Diagonal(x=expDt$UNC^2),S_pwl,P_pwl)
updSysDt[ERRTYPE=="pw", V1 := as.vector(polygon_chain_parameters)]

# posterior mean of the polygon chain parameters given the experimental data
# including the MLO added uncertainties
polygon_chain_parameters <- prior_polygon_chain_parameters +
														t(S_pwl %*% P_pwl) %*% mult_invCov_x(expDt[,DATA]-prior_exp_data,Diagonal(x=expDt$UNC^2),S_extra,U_extra)

updSysDt[ERRTYPE=="pw", V2 := as.vector(polygon_chain_parameters)]


for (curReac in unique(updSysDt[ERRTYPE=="pw"]$EXPID)) {

	curSysDt <- updSysDt[ERRTYPE=="pw" & EXPID==curReac]
	v <- curSysDt$DATA

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

	#the prior
	v <- sysDt[ERRTYPE=="pw" & EXPID==curReac]$DATA
	y <- Z %*% v
	updSysDt[ERRTYPE=="pw" & EXPID==curReac,XS_POLY_CHAIN_PRIOR:=as.vector(y)]
}

updSysDt[ERRTYPE=="pw",EN_POLY_CHAIN:=EnGridPolyChains]

#updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,INL-01",REAC:="(24-CR-52(N,INL)24-CR-52,,SIG)"]
#updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,P-01",REAC:="(24-CR-52(N,P)23-V-52,,SIG)"]
#updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,2N-01",REAC:="(24-CR-52(N,2N)24-CR-51,,SIG)"]
#updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,TOT-01",REAC:="(24-CR-52(N,TOT),,SIG)"]
#updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,EL-01",REAC:="(24-CR-52(N,EL)24-CR-52,,SIG)"]
#updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,A-01",REAC:="(24-CR-52(N,A)22-TI-49,,SIG)"]
#updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,D-01",REAC:="(24-CR-52(N,D)23-V-51,,SIG)"]

#updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,X-01",REAC:="(24-CR-52(N,X)1-H-1,,SIG)"]
#updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,X-02",REAC:="(24-CR-52(N,X)2-HE-4,,SIG)"]
mapAssignment <- reacHandler$getMapAssignment()
for(i in seq_len(nrow(mapAssignment))) {
  updSysDt[ERRTYPE=="pw" & EXPID==mapAssignment$EXPID[i],
           REAC := mapAssignment$REAC[i]]
}

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
#geom_line(data = updSysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN_PRIOR, col='prior'), linewidth = 0.1) +
scale_color_manual(name='Regression Model',
                     breaks=c('orig. sys. unc.', 'only stat. unc.', 'extra sys. unc.','prior'),
                     values=c('extra sys. unc.'='green', 'only stat. unc.'='red', 'orig. sys. unc.'='black','prior'='blue')) +
ylim(0,NA) +
facet_wrap(~ REAC, ncol=3, scales='free')

dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
filepath <- file.path(plotPath, 'MLO_uncertainties.png')
#ggsave(filepath, ggp1, width = 53, height = 29.8125, units = "cm", dpi = 300)
ggsave(filepath, ggp1, width = 0.2*53, height = 0.2*29.8125, units = "cm", dpi = 300)

plotDir <- file.path(plotPath, 'MLO_uncertainties')
dir.create(plotDir, recursive=TRUE, showWarnings=FALSE)
for(channel in unique(expDt[,REAC])) {
	cur_expDt <- expDt[REAC==channel]
	cur_SysDt <- updSysDt[REAC==channel]

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
	#geom_line(data = cur_SysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN_PRIOR, col='prior'), linewidth = 0.1) +
	scale_color_manual(name='Regression Model',
	                     breaks=c('orig. sys. unc.', 'only stat. unc.', 'extra sys. unc.','prior'),
	                     values=c('extra sys. unc.'='green', 'only stat. unc.'='red', 'orig. sys. unc.'='black','prior'='blue')) +
	ylim(0,NA)

	filename <- file.path(plotDir, paste0(channel,'.png'))
	ggsave(filename, ggp1, width = 0.2*53, height = 0.2*29.8125, units = "cm", dpi = 300)
}
