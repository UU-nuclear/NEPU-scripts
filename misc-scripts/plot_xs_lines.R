library(ggplot2)
#library(mvtnorm)
library(mvnfast)

needsDt <- read_object(1, "needsDt")
expDt <- read_object(3, "expDt")

finalPars <- read_object(11, "finalPars")
finalParCovmat <- read_object(11, "finalParCovmat")
allResults <- read_object(12, "allResults")
allParsets <- read_object(12, "allParsets")

reactions <- expDt[,unique(REAC)]

# create a data.table with the model energy-grid
# and the reaction channel specification

extNeedsDt <- needsDt[,{
    stopifnot(all(L2 == 0) & all(L3 == 0))
    list(L1 = defineEnergyGrid(L1, energyGridrandomFiles, enPolicy="compgrid"),
         L2 = 0, L3 = 0)
}, by=c("PROJECTILE", "ELEMENT", "MASS", "REAC")]
extNeedsDt[, IDX := seq_len(.N)]

# create the grid which to interpolate to

# create the interpolation matrix using the exforHandler
modSubents <- lapply(reactions, createSubentStub, en=energyGridrandomFiles) # fake exfor sub-entries
modDt <- exforHandler$extractData(modSubents, ret.values=FALSE) # fake experimental data
Smod <- exforHandler$getJac(modDt, extNeedsDt, modSubents)
setkey(modDt, IDX, DIDX)

modDt[,L2:=NULL]
modDt[,L3:=NULL]

#####################################
# Get the estimated posterior probability according to the MVN distribution
# I only consider the probability for the fitted parameters, for now.
# I am not considering the Gaussian Process extrapolated parameters

post_probs <- dmvn(t(allParsets[1:length(finalPars),]),
					mu=finalPars,
					sigma = finalParCovmat,
					log = TRUE,
					ncores = 1,
					isChol = FALSE)

# these are needed to calculate the exact posterior prob.
optParamDt <- read_object(10, "optParamDt")
# P0 <- read_object(10, "P0") This one is only for the adjustable parameters
optSysDt <- read_object(6, "optSysDt")
optGpDt <- read_object(6, "optGpDt")
refParamDt <- read_object(2, "refParamDt")
extNeedsDt2 <- read_object(2, "extNeedsDt")
modList <- read_object(3, "modList")
fullSensDt <- read_object(5, "fullSensDt") 
subents <- read_object(1, "subents")

# reconstruct the covariance matrix P containing both the 
# covariances of the systematic experimental errors and
# the covariances for the model parameters. This matrix is diagonal
# for the portions relating to experiments and energy-independent
# parameters. The blocks associated with energy-dependent TALYS
# parameters are constructed by a Gaussian process and contain
# therefore correlations.
# prepare the Gaussian process handler 
gpHandler <- createSysCompGPHandler()

# prepare the handlers to map systematic uncertainties of the experiments
normHandler <- createSysCompNormHandler("DATAREF")
normHandler$addSysUnc("EXPID", "", 0, 0, FALSE)

# prepare the TALYS handler to map from model parameters to predictions
talysHandler <- createSysCompModelHandler()
talysHandler$setRef(extNeedsDt2, fullSensDt, refParamDt,
                    exforHandler, c(subents, modList$SUBENT))
talysHandler$setPrior(refParamDt)

# create global handler and register the individual handlers
sysCompHandler <- createSysCompHandler()
sysCompHandler$addHandler(normHandler)
sysCompHandler$addHandler(talysHandler)
sysCompHandler$addGPHandler(gpHandler)
P <- sysCompHandler$cov(optSysDt, optGpDt, ret.mat = TRUE)
expSel <- optSysDt[, !grepl("TALYS-", EXPID)]
talysSel <- !expSel

# Prior parameter covariance matrix of all parameters
P0 <- P[talysSel, talysSel]

# sensitivty matrix to interpolate the calculations
Sexp <- exforHandler$getJac(expDt, extNeedsDt, subents)


post_probs2 <- rep(0,ncol(allParsets))
inv_finalParCovmat <- solve(finalParCovmat)
tmp <- determinant(finalParCovmat)
stopifnot(tmp$sign == 1)
logDetfinalParCovmat <- as.vector(tmp$modulus)

for(calc in seq_len(ncol(allParsets)))
{
	# interpolate talys to the experimental energies
	f_sample <- Sexp %*% allResults[,calc]

	dfinal <- allParsets[1:length(finalPars),calc] - finalPars

	chisqr <- as.vector( t(dfinal) %*%  inv_finalParCovmat %*% dfinal)

  	post_probs2[calc] <- - 0.5 * (
  		chisqr + logDetfinalParCovmat + 
  		length(finalPars)* log(2*pi)
  	)
}

# post_probs2 == post_probs should hold
cat("max(abs(post_probs - post_probs2)) = ", max(abs(post_probs - post_probs2)), "\n")
#####################################
# Calculate the actual posterior probability for each sample
# I only consider the probability for the fitted parameters, for now.
# I am not considering the Gaussian Process extrapolated parameters

# experimental quantities
yexp <- read_object(10, "yexp")
D <- read_object(10, "D")
S <- read_object(10, "S0")
X <- read_object(10, "X")
cholZ <- makeCholZ(D, S, X)

# prior parameter mean
priorPar <- unlist(optParamDt[3:nrow(optParamDt),]$PARVAL) 

# inverse of prior covariance matrix
invP0 <- solve(P0) 

# this should be done but is not in the fit (would make no difference if all == 1)
# priorPar <- paramTrafo$fun(priorPar) 

post_probs_real <- rep(0,ncol(allParsets))

logDetExp <- logDetCov(D, S, P, cholZ)

tmp <- determinant(P0)
stopifnot(tmp$sign == 1)
logDetPar <- as.vector(tmp$modulus)

for(calc in seq_len(ncol(allParsets)))
{
	# interpolate talys to the experimental energies
	f_sample <- Sexp %*% allResults[,calc]

	dprior <- allParsets[1:length(priorPar),calc] - priorPar

	Lprior <- as.vector(crossprod(dprior, invP0 %*% dprior))
	#cat(Lprior - as.vector(t(dprior) %*% invP0 %*% dprior),"\n")
	chisqr <- as.vector(mult_xt_invCov_x(yexp - f_sample, D, S, X, cholZ = cholZ))
	#post_probs_real[calc] <- -0.5 * (as.vector(mult_xt_invCov_x(yexp - f_sample, D, S, X, cholZ = cholZ)) + Lprior)

  	post_probs_real[calc] <- - 0.5 * (
  		chisqr + logDetExp + 
  		Lprior + logDetPar + 
  		(length(yexp) + length(priorPar))* log(2*pi)
  	)
}
#####################################

defaultMinEnergy <- 1
maxEnergy <- 200

# reoder the data so that the REAC keyword is correct
reorder <- function(x)
{
	as.vector(Smod %*% x)
}
allResults <- apply(allResults,2,reorder)

samplesDt <- data.table(L1=modDt$L1,
						REAC=modDt$REAC,
						DATA = as.vector(allResults[,2:ncol(allResults)]),
						#LOGP = rep(post_probs[2:length(post_probs)],each=length(modDt$L1)),
						LOGP = rep(post_probs2[2:length(post_probs2)] - post_probs2[1],each=length(modDt$L1)),
						LOGP_REAL = rep(post_probs_real[2:length(post_probs_real)]-post_probs_real[1],each=length(modDt$L1))
						#LOGP_REAL = rep(post_probs_real[2:length(post_probs_real)],each=length(modDt$L1))
						)

modeDt <- data.table(L1=modDt$L1,
						REAC=modDt$REAC,
						DATA = as.vector(allResults[,1]),
						#LOGP = rep(post_probs[1],length(modDt$L1)),
						LOGP = rep(post_probs2[1],length(modDt$L1)),
						LOGP_REAL = rep(post_probs_real[1],length(modDt$L1))
						)

plot_MVNprob <- ggplot(data=samplesDt[L1<maxEnergy],mapping = aes(x=L1,y=DATA)) + theme_bw() +
		geom_line(aes(col=LOGP,group=LOGP)) +
		geom_line(data=modeDt[L1<maxEnergy],col='red') +
		facet_wrap(~REAC,scales="free_y")
print(plot_MVNprob)

plot_REALprob <- ggplot(data=samplesDt[L1<maxEnergy],mapping = aes(x=L1,y=DATA)) + theme_bw() +
		geom_line(aes(col=LOGP_REAL,group=LOGP_REAL)) +
		geom_line(data=modeDt[L1<maxEnergy],col='red') +
		#geom_point(data=expDt[L1<maxEnergy],alpha=0.25) +
		#geom_errorbar(data=expDt[L1<maxEnergy],aes(ymin=DATA-UNC,ymax=DATA+UNC),alpha=0.25) +
		facet_wrap(~REAC,scales="free_y")
print(plot_REALprob)

# This should tell us something about how good our approximation is
# For a perfect approximation this ratio is always 1
hist((post_probs-post_probs[1]) / (post_probs_real-post_probs_real[1]),breaks=100)
#hist((post_probs-post_probs[1]) - (post_probs_real-post_probs_real[1]),breaks=100)

# log(P(p)) - log(P(p_max)) : pdf normalized to 1 at the maximum 

#dmvn(t(allParsets[1:length(finalPars),]),

#test <- as.vector(finalPars)
#test <- test + runif(length(test),-0.001,0.001)
#dmvn(test,
#					mu=finalPars,
#					sigma = finalParCovmat,
#					log = TRUE,
#					ncores = 1,
#					isChol = FALSE)