#################################################
#       SCRIPT Setup
##################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}

library(ggplot2)
library(latex2exp)

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

subents <- read_object(1, "subents")
needsDt <- read_object(1, "needsDt")
extNeedsDt <- read_object(2, "extNeedsDt")
origSysDt <- read_object(4, "origSysDt")
updSysDt <- read_object(4, "updSysDt")
modDt <- read_object(3,"modDt")
optExpDt <- read_object(7, "optExpDt")
allResults <- read_object(12, 'allResults')

optParamDt <- read_object(10,"optParamDt")
optRes <- read_object(10, "optRes")
P0 <- read_object(10, "P0")
X <- read_object(10, "X")
SX <- read_object(10, "S0")

# sort the experimental data table by recation, energy
optExpDt <- optExpDt[order(REAC,L1)]
optExpDt[, IDX := seq_len(.N)]

reactions <- optExpDt[,unique(REAC)]

normHandler <- createSysCompNormHandler("DATAREF")
normHandler$addSysUnc("EXPID", "", 0, 0, TRUE)

sysCompHandler <- createSysCompHandler()
sysCompHandler$addHandler(normHandler)

S <- sysCompHandler$map(optExpDt, origSysDt, ret.mat = TRUE)
origX <- sysCompHandler$cov(origSysDt, ret.mat = TRUE)
updX <- sysCompHandler$cov(updSysDt, ret.mat = TRUE)
statUnc <- getDt_UNC(optExpDt)

origUnc <- sqrt(statUnc^2 + diag(S %*% origX %*% t(S))) 
updUnc <- sqrt(statUnc^2 + diag(S %*% updX %*% t(S)))

setkey(optExpDt, IDX)
optExpDt[, ORIGUNC := origUnc]
optExpDt[, UPDUNC := updUnc]

# ----------Calculate the Chi-squares -----------
RandomFilesNeedsDt <- needsDt[,{
    stopifnot(all(L2 == 0) & all(L3 == 0))
    list(L1 = defineEnergyGrid(L1, energyGridrandomFiles, enPolicy="compgrid"),
         L2 = 0, L3 = 0)
}, by=c("PROJECTILE", "ELEMENT", "MASS", "REAC")]
RandomFilesNeedsDt[, IDX := seq_len(.N)]

allResults[is.na(allResults)] <- 0
sampledResults <- allResults[,2:ncol(allResults)]

Sexp <- exforHandler$getJac(optExpDt,RandomFilesNeedsDt,subents)
transform_function <- function(x) { as.vector(Sexp%*%x) }
allResults_transform <- apply(allResults,2,FUN=transform_function)
sampledResults_transform <- allResults_transform[,2:ncol(allResults_transform)]

uncinfo_transform <- cov.wt(t(sampledResults_transform))
model_covariance <- uncinfo_transform$cov
model_predictions <- uncinfo_transform$center
model_uncertainties <- sqrt(diag(model_covariance))
optExpDt <- optExpDt[order(REAC,L1)]
optExpDt[, FIT_MEAN := model_predictions]

optExpDt[, FIT_MODE := allResults_transform[,1]]


stat_unc <- getDt_UNC(optExpDt)
D <- Diagonal(x = stat_unc^2)

S <- sysCompHandler$map(optExpDt, origSysDt, ret.mat = TRUE)
P <- sysCompHandler$cov(updSysDt, ret.mat = TRUE)
d <- optExpDt[,DATA-FIT_MEAN]
chi2 <- chisquare(d,D,S,P)

n_points <- nrow(optExpDt)
n_pars_energy_indep <- nrow(optParamDt[ADJUSTABLE==TRUE & !grepl("\\(.+\\)",PARNAME)])
n_pars_energy_dep <- length(optParamDt[ADJUSTABLE==TRUE & grepl("\\(.+\\)",PARNAME),unique(str_remove(PARNAME,"\\(.+\\)"))])

# I consider two quantities derived from the data for each energy dependent parameter (length-scale, amplitude)-hyper-parameters

degrees_of_freedom <- n_points - (n_pars_energy_indep + 2*n_pars_energy_dep)

cat("chi-square/degree-of-freedom = ",chi2," / ",degrees_of_freedom," = ",chi2/degrees_of_freedom,"\n")

cat("chi-square p-value = ",pchisq(chi2,degrees_of_freedom),"\n")

# calculate also the Chi2 based on the derived model covariance
# I can't do it this way, because of the linear interpolation
# the matrix becomes exactly singular
# chi2 <- t(d) %*% solve(model_covariance,d)
# chi2 <- t(d) %*% Matrix::solve(model_covariance,d)