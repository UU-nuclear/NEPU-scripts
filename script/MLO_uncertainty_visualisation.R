library(ggnewscale)

subents <- read_object(1, "subents")
modList <- read_object(3, "modList")
modDt <- read_object(3, "modDt") # default model prediction mapped to experimental energies
expDt <- read_object(3, "expDt")
sysUncDt <- read_object(3, "sysUncDt")

#subents <- read_object(1, "orig_subents")
#modList <- read_object(3, "orig_modList")
#modDt <- read_object(3, "orig_modDt") # default model prediction mapped to experimental energies
#expDt <- read_object(3, "orig_expDt")
#sysUncDt <- read_object(3, "sysUncDt")

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

# the maximum energy is taken from the limits defined in the config file
#En_max <- maxExpEn + 2

curReac <- "(26-FE-56(N,INL)26-FE-56,,SIG)"
En_max <- max(expDt[REAC==curReac]$L1) + 0.1
curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), En_max, by = 0.1)
curUncs <- c(1e4, 1e4, rep(2, length(curEnGrid)-2))
reacHandler$assignMapToReac("pw", curReac,
                            vals = rep(0, length(curEnGrid)),
                            uncs = curUncs,
                            opts = list(ens = curEnGrid,
                                        order = 3, outsideZero = TRUE))

tot_length <- length(curEnGrid)
EnGridPolyChains <- curEnGrid

curReac <- "(26-FE-56(N,P)25-MN-56,,SIG)"
En_max <- max(expDt[REAC==curReac]$L1) + 0.1
curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), En_max, by = 0.1)
curUncs <- c(1e4, 1e4, rep(2, length(curEnGrid)-2))
reacHandler$assignMapToReac("pw", curReac,
                            vals = rep(0, length(curEnGrid)),
                            uncs = curUncs,
                            opts = list(ens = curEnGrid,
                                        order = 3, outsideZero = TRUE))

tot_length <- tot_length + length(curEnGrid)
EnGridPolyChains <- c(EnGridPolyChains,curEnGrid)

curReac <- "(26-FE-56(N,2N)26-FE-55,,SIG)"
En_max <- max(expDt[REAC==curReac]$L1) + 0.1
curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), En_max, by = 0.1)
curUncs <- c(1e4, 1e4, rep(0.5, length(curEnGrid)-2))
reacHandler$assignMapToReac("pw", curReac,
                            vals = rep(0, length(curEnGrid)),
                            uncs = curUncs,
                            opts = list(ens = curEnGrid,
                                        order = 3, outsideZero = TRUE))

tot_length <- tot_length + length(curEnGrid)
EnGridPolyChains <- c(EnGridPolyChains,curEnGrid)

curReac <- "(26-FE-56(N,TOT),,SIG)"
bin_width <- 0.3
En_max <- max(expDt[REAC==curReac]$L1) + bin_width
curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), En_max, by = bin_width)
#curEnGrid <- seq(0, En_max, by = 0.1)
#En_min <- min(expDt[REAC==curReac]$L1) - 0.1
#curEnGrid <- seq(En_min, En_max, by = 0.1)
curUncs <- c(1e5, 1e5, rep(20, length(curEnGrid)-2))
reacHandler$assignMapToReac("pw", curReac,
                            #vals = rep(0, length(curEnGrid)),
                            vals = c(0,rep(0, length(curEnGrid)-1)),
                            uncs = curUncs,
                            opts = list(ens = curEnGrid,
                                        order = 3, outsideZero = TRUE))

tot_length <- tot_length + length(curEnGrid)
EnGridPolyChains <- c(EnGridPolyChains,curEnGrid)


curReac <- "(26-FE-56(N,EL)26-FE-56,,SIG)"
En_max <- max(expDt[REAC==curReac]$L1) + 0.1
curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), En_max, by = 0.1)
curUncs <- c(1e4, 1e4, rep(2, length(curEnGrid)-2))
reacHandler$assignMapToReac("pw", curReac,
                            vals = rep(0, length(curEnGrid)),
                            uncs = curUncs,
                            opts = list(ens = curEnGrid,
                                        order = 3, outsideZero = TRUE))

tot_length <- tot_length + length(curEnGrid)
EnGridPolyChains <- c(EnGridPolyChains,curEnGrid)


curReac <- "(26-FE-56(N,A)24-CR-53,,SIG)"
En_max <- max(expDt[REAC==curReac]$L1) + 0.1
curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), En_max, by = 0.1)
curUncs <- c(1e4, 1e4, rep(50, length(curEnGrid)-2))
reacHandler$assignMapToReac("pw", curReac,
                            vals = rep(0, length(curEnGrid)),
                            uncs = curUncs, 
                            opts = list(ens = curEnGrid,
                                        order = 3, outsideZero = TRUE))

tot_length <- tot_length + length(curEnGrid)
EnGridPolyChains <- c(EnGridPolyChains,curEnGrid)


curReac <- "(26-FE-56(N,D)25-MN-55,,SIG)"
En_max <- max(expDt[REAC==curReac]$L1) + 0.1
curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), En_max, by = 0.1)
curUncs <- c(1e4, 1e4, rep(5, length(curEnGrid)-2))
reacHandler$assignMapToReac("pw", curReac,
                            vals = rep(0, length(curEnGrid)),
                            uncs = curUncs, 
                            opts = list(ens = curEnGrid,
                                        order = 3, outsideZero = TRUE))

tot_length <- tot_length + length(curEnGrid)
EnGridPolyChains <- c(EnGridPolyChains,curEnGrid)

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

# Fix some experimental normalization uncertainties 
# to help MLO find a good solution

# (26-FE-56(N,INL)26-FE-56,,SIG)"
curExpId <- "23134005"
#sysDt[paste0("EXPID-",curExpId) == EXPID, UNC := 0.005]
#expDt[curExpId == EXPID, UNC := pmax(DATAREF * 0.05, 1)]
#sysDt[paste0("EXPID-",curExpId) == EXPID, ADJUSTABLE := FALSE]

# ((26-FE-56(N,2N)26-FE-55,,SIG)
curExpId <- "23171003"
#sysDt[paste0("EXPID-",curExpId) == EXPID, UNC := 0.005]
#expDt[curExpId == EXPID, UNC := pmax(DATAREF * 0.05, 1)]
#sysDt[paste0("EXPID-",curExpId) == EXPID, ADJUSTABLE := FALSE]

# ----------------------------------------------------------------------------------
# Condition the polygon-chain GP on the experiemtnal data
# and set the DATA column in sysDt

## mapping matrix for linear interpolation of the polygon chain
## the UNC column is the same in the original and the updated
## SysDt[ERRTYPE=="pw"]
#S_pwl <- sysCompHandler$map(expDt, sysDt[ERRTYPE=="pw"], ret.mat = TRUE)
## prior on the second derivative of the polygon chain
#P_pwl <- sysCompHandler$cov(sysDt[ERRTYPE=="pw"], ret.mat = TRUE)
#
## mapping matrix for systematics to experimental data
#S <- sysCompHandler$map(expDt, sysDt, ret.mat = TRUE)
## Systematic uncertainties of experimental data (before MLO) + 
## prior on the second derivative of the polygon chain
#U <- sysCompHandler$cov(sysDt, ret.mat = TRUE)
#
## posterior mean of the polygon chain parameters given the experimental data
## and their covariance
#polygon_chain_parameters <- t(S_pwl %*% P_pwl) %*% mult_invCov_x(getDt_DATA(expDt),Diagonal(x=expDt$UNC^2),S,U)
##polygon_chain_parameters <- t(S %*% U) %*% mult_invCov_x(getDt_DATA(expDt),Diagonal(x=expDt$UNC^2),S,U)
#sysDt[ERRTYPE=="pw", DATA := as.vector(polygon_chain_parameters)]
# ----------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------
# Condition the polygon-chain GP on the experimental data
# and set the DATA column in sysDt (using only the stat. uncs)

# mapping matrix for linear interpolation of the polygon chain
# the UNC column is the same in the original and the updated
# SysDt[ERRTYPE=="pw"]
S_pwl <- sysCompHandler$map(expDt, sysDt[ERRTYPE=="pw"], ret.mat = TRUE)
# prior on the second derivative of the polygon chain
P_pwl <- sysCompHandler$cov(sysDt[ERRTYPE=="pw"], ret.mat = TRUE)
# posterior mean of the polygon chain parameters given the experimental data
# ignoring the systematic errors in the data
polygon_chain_parameters <- t(S_pwl %*% P_pwl) %*% mult_invCov_x(getDt_DATA(expDt),Diagonal(x=expDt$UNC^2),S_pwl,P_pwl)
sysDt[ERRTYPE=="pw", DATA := as.vector(polygon_chain_parameters)]
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
	#initUncs <- lowerLims + runif(length(lowerLims)) * (upperLims-lowerLims)
	#initUncs <- lowerLims
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

# posterior mean of the polygon chain parameters given the experimental data
# and their covariance
polygon_chain_parameters <- t(S_pwl %*% P_pwl) %*% mult_invCov_x(getDt_DATA(expDt),Diagonal(x=expDt$UNC^2),S,U)
#polygon_chain_parameters <- t(S %*% U) %*% mult_invCov_x(getDt_DATA(expDt),Diagonal(x=expDt$UNC^2),S,U)
updSysDt[ERRTYPE=="pw", DATA := as.vector(polygon_chain_parameters)]

# posterior mean of the polygon chain parameters given the experimental data
# ignoring the systematic errors in the data
polygon_chain_parameters <- t(S_pwl %*% P_pwl) %*% mult_invCov_x(getDt_DATA(expDt),Diagonal(x=expDt$UNC^2),S_pwl,P_pwl)
updSysDt[ERRTYPE=="pw", V1 := as.vector(polygon_chain_parameters)]

# posterior mean of the polygon chain parameters given the experimental data
# including the MLO added uncertainties
polygon_chain_parameters <- t(S_pwl %*% P_pwl) %*% mult_invCov_x(getDt_DATA(expDt),Diagonal(x=expDt$UNC^2),S_extra,U_extra)
#polygon_chain_parameters <- t(S_extra %*% U_extra) %*% mult_invCov_x(getDt_DATA(expDt),Diagonal(x=expDt$UNC^2),S_extra,U_extra)
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
}

updSysDt[ERRTYPE=="pw",EN_POLY_CHAIN:=EnGridPolyChains]

updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,INL-01",REAC:="(26-FE-56(N,INL)26-FE-56,,SIG)"]
updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,P-01",REAC:="(26-FE-56(N,P)25-MN-56,,SIG)"]
updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,2N-01",REAC:="(26-FE-56(N,2N)26-FE-55,,SIG)"]
updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,TOT-01",REAC:="(26-FE-56(N,TOT),,SIG)"]
updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,EL-01",REAC:="(26-FE-56(N,EL)26-FE-56,,SIG)"]
updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,A-01",REAC:="(26-FE-56(N,A)24-CR-53,,SIG)"]
updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,D-01",REAC:="(26-FE-56(N,D)25-MN-55,,SIG)"]

ggp1 <- ggplot(expDt) + theme_bw() +
theme(axis.text=element_text(size=9),
	               axis.title=element_text(size=10),
	               plot.title=element_text(size=12)) +
xlab("energy [MeV]") + ylab("cross section [mbarn]") +
guides(col = "none") +
geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = "black",
                               size = 0.5, width = 0.3) +
geom_errorbar(aes(x = L1, ymin = DATA - ORIGUNC, ymax = DATA + ORIGUNC, col = EXPID),
                           size = 0.5, width = 0.2) +
geom_point(aes(x = L1, y = DATA, col = EXPID)) +

new_scale_colour() +

geom_line(data = updSysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN, col='orig. sys. unc.')) +
geom_line(data = updSysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN_STAT, col='only stat. unc.')) +
geom_line(data = updSysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN_MLO, col='extra sys. unc.')) +
scale_color_manual(name='Regression Model',
                     breaks=c('orig. sys. unc.', 'only stat. unc.', 'extra sys. unc.'),
                     values=c('extra sys. unc.'='green', 'only stat. unc.'='red', 'orig. sys. unc.'='black')) +
facet_wrap(~ REAC, ncol=3, scales='free')

#ggsave(file.path(plotPath, 'MLO.png'), ggp1)
#dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
#filepath <- file.path(plotPath, 'piece_wise_linear_fit.png')
#ggsave(filepath, ggp, width = 53, height = 29.8125, units = "cm", dpi = 300)

#ggp <- ggplot(expDt[REAC=="(26-FE-56(N,TOT),,SIG)" & EXPID!=22316003 & EXPID!=13764002]) + theme_bw() +
#ggp <- ggplot(expDt[REAC=="(26-FE-56(N,TOT),,SIG)" & EXPID!=22316003]) + theme_bw() +
#ggp <- ggplot(expDt[REAC=="(26-FE-56(N,TOT),,SIG)" & EXPID!=13764002]) + theme_bw() +
#ggp2 <- ggplot(expDt[REAC=="(26-FE-56(N,TOT),,SIG)"]) + theme_bw() +
ggp2 <- ggplot(expDt[REAC=="(26-FE-56(N,2N)26-FE-55,,SIG)"]) + theme_bw() +
theme(axis.text=element_text(size=9),
	               axis.title=element_text(size=10),
	               plot.title=element_text(size=12)) +
xlab("energy [MeV]") + ylab("cross section [mbarn]") +
#geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = "black",
#                               size = 0.5, width = 0.3) +
geom_errorbar(aes(x = L1, ymin = DATA - ORIGUNC, ymax = DATA + ORIGUNC, col = EXPID),
                           size = 0.5, width = 0.2) +
geom_point(aes(x = L1, y = DATA, col = EXPID)) +

new_scale_colour() +
#geom_line(data = updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,TOT-01"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN, col='orig. sys. unc.')) +
geom_line(data = updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,2N-01"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN_STAT, col='only stat. unc.')) +
#geom_line(data = updSysDt[ERRTYPE=="pw" & EXPID=="REACEXP-N,TOT-01"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN_MLO, col='extra sys. unc.')) +
scale_color_manual(name='Regression Model',
                     breaks=c('orig. sys. unc.', 'only stat. unc.', 'extra sys. unc.'),
                     values=c('extra sys. unc.'='green', 'only stat. unc.'='red', 'orig. sys. unc.'='black')) +
geom_line(aes(x = L1, y = DATAREF),col="blue")

#ggsave(file.path(plotPath, 'MLO_unc_vis.png'), ggp1)

#Stest <- sysCompHandler2$map(expDt[REAC=="(26-FE-56(N,TOT),,SIG)"], origSysDt, ret.mat = TRUE)
#origX <- sysCompHandler2$cov(origSysDt, ret.mat = TRUE)
#updX <- sysCompHandler2$cov(updSysDt, ret.mat = TRUE)
#
#library(reshape2)
#library(gplots)
#library(RColorBrewer)
#
#covMat1 <- S %*% origX %*% t(S)
#color_palette <- colorRampPalette(brewer.pal(8, "Blues"))(25)
#heatmap.2(as.matrix(covMat1),Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,
#	scale="none",margins=c(8,8),trace = "none",
#	dendrogram="none",density.info = "none")

# If I do the following:
# P <- sysCompHandler$cov(origSysDt, ret.mat = TRUE)
# S <- sysCompHandler$map(expDt, origSysDt, ret.mat = TRUE)
#
# S_pwl <- sysCompHandler$map(expDt, origSysDt[ERRTYPE=="pw"], ret.mat = TRUE)
# P_pwl <- sysCompHandler$cov(origSysDt[ERRTYPE=="pw"], ret.mat = TRUE)
#
# v2 <- t(S_pwl %*% P_pwl) %*% mult_invCov_x(getDt_DATA(expDt),Diagonal(x=expDt$UNC^2),S,U)
# v1 <- t(S %*% P) %*% mult_invCov_x(getDt_DATA(expDt),Diagonal(x=expDt$UNC^2),S,U)
#
# then 
#
# all(v1[1:length(v2)] == v2) is true
#
# That is, I get the same posterior for the polygon chain parameters from the two v1 and v2 above
#
