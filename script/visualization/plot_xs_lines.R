library(ggplot2)
#library(mvtnorm)
library(mvnfast)
library(latex2exp)

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

#####################################
# Get the estimated posterior probability according to the MVN distribution,
# but using the GLS approximation (instead of Georgs)
# I only consider the probability for the fitted parameters, for now.
# I am not considering the Gaussian Process extrapolated parameters
optRes <- read_object(10, "optRes")
D <- read_object(10, "D")
S0 <- read_object(10, "S0")
X <- read_object(10, "X")
optSysDt_optpars <- read_object(10, "optSysDt_optpars")
optSysDt_allpars <- read_object(10, "optSysDt_allpars")
optGpDt <- read_object(6, "optGpDt")

# compute the GLS Hessian
# Smod maps from model parameters to optimzed parameters only
P0_all <- sysCompHandler$cov(optSysDt_allpars, optGpDt, ret.mat = TRUE)

Smod <- optRes$jac
tS_invCexp_S <- mult_xt_invCov_x(Smod, D, S0, X)
invP0_all <- solve(P0_all)
H_gls <- (-invP0_all)
optpars_indices <- optSysDt_optpars[, sort(IDX)]
H_gls[optpars_indices, optpars_indices] <- H_gls[optpars_indices, optpars_indices] - (tS_invCexp_S)
Covmat_GLS <- (-1) * solve(H_gls)

post_probs_GLS <- dmvn(t(allParsets[1:length(finalPars),]),
					mu=finalPars,
					sigma = Covmat_GLS,
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
priorPar <- unlist(refParamDt[ADJUSTABLE==TRUE,PARVAL]) 

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
		geom_line(aes(col=LOGP,group=LOGP),size=0.25) +
		geom_point(data=expDt[L1<maxEnergy],alpha=0.25,size=0.1) +
		geom_errorbar(data=expDt[L1<maxEnergy],aes(ymin=DATA-UNC,ymax=DATA+UNC),alpha=0.25) +
		geom_line(data=modeDt[L1<maxEnergy],col='red',size=0.25) +
		facet_wrap(~REAC,scales="free_y")
print(plot_MVNprob)

plot_REALprob <- ggplot(data=samplesDt[L1<maxEnergy],mapping = aes(x=L1,y=DATA)) + theme_bw() +
		geom_line(aes(col=LOGP_REAL,group=LOGP_REAL),size=0.25) +
		geom_point(data=expDt[L1<maxEnergy],alpha=0.25,size=0.1) +
		geom_errorbar(data=expDt[L1<maxEnergy],aes(ymin=DATA-UNC,ymax=DATA+UNC),alpha=0.25) +
		geom_line(data=modeDt[L1<maxEnergy],col='red',size=0.25) +
		facet_wrap(~REAC,scales="free_y")
print(plot_REALprob)

dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
filepath <- file.path(plotPath, paste0('sampled_xs_MVNprob.png'))
ggsave(filepath, plot_MVNprob, width = 16*1.5, height = 9*1.5, units = "cm", dpi = 300)
filepath <- file.path(plotPath, paste0('sampled_xs_MVNprob2.png'))
ggsave(filepath, plot_MVNprob+xlim(0,5), width = 16*1.5, height = 9*1.5, units = "cm", dpi = 300)

dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
filepath <- file.path(plotPath, paste0('sampled_xs_REALprob.png'))
ggsave(filepath, plot_REALprob, width = 16*1.5, height = 9*1.5, units = "cm", dpi = 300)
filepath <- file.path(plotPath, paste0('sampled_xs_REALprob2.png'))
ggsave(filepath, plot_REALprob+xlim(0,5), width = 16*1.5, height = 9*1.5, units = "cm", dpi = 300)

# check how good our approximated posterior pdf is:

# We don't know the integral of the real pdf, so we can't compare the real probability
# with the MVN approxiamtion. However, we do know the maximum probability density of
# both the real and the MVN approximation (the first 'sample'), so we can make a realtive comparison.
#    let P_k be the probability of a sample k
#    let ln(P_k) be the logarithm of the probability
# we know that  post_probs_real[k] = Norm*ln(P_k), where norm is a normalization constant, so
#    yy = post_probs_real[k] - post_probs_real[1] = ln(P_k/P_1)
# the same can be done for the MVN approximation Pi of the prob.
#    xx = post_probs[k] - post_probs[1] = ln(Pi_k/Pi_1)


# The difference xx - yy = ln((P_k/P_1)/(Pi_k/Pi_1))



probDt <- data.table(xx = post_probs[2:length(post_probs)],
	yy = post_probs_real[2:length(post_probs_real)])
ggplot(data=probDt, mapping=aes(x=xx,y=yy)) + theme_bw() +
	geom_point() +
	ylab(bquote(ln(p[true])+ln(N))) +
	xlab(bquote(ln(p[approx]))) + 
	geom_smooth(method='lm', formula= y~x)

# The following presentation of the data is pretty good, it shows the renormalized
# (to the height of the pdfs) logarithm of the approximated vs. the true prob. density.
# Perfect approximation would mean all points lie on the diagonal, points above the diagonal
# means the true probability is larger than the estimated one, while points below the line
# means the true probability is smaller than the estimated one.
# If we see more points above the line our estimate is not conservative enough
# If we see more points below the line our estimate is (too) conservative

yy <- post_probs_real[2:length(post_probs_real)] - post_probs_real[1]
xx <- post_probs[2:length(post_probs)] - post_probs[1]
zz <- post_probs_GLS[2:length(post_probs_GLS)] - post_probs_GLS[1]

probDt2 <- data.table(xx = xx,
	yy = yy, zz = zz)
ggplot(data=probDt2, mapping=aes(x=xx,y=yy)) + theme_bw() +
	geom_point() +
	geom_point(mapping=aes(x=zz,y=yy),col='red') +
	xlab(TeX("$\\propto \\ln(P_{approx.})$")) +
	ylab(TeX("$\\propto \\ln(P_{true})$")) +
	geom_abline()

# I should add the GLS estimate of the probability density!!!!