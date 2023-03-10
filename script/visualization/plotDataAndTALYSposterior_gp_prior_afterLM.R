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

getChi2reaction <- function(reaction_string) {
  reaction.expDt <- copy(optExpDt[REAC==reaction_string])
  reaction.expDt[, IDX := seq_len(.N)]
  reaction.stat_unc <- getDt_UNC(reaction.expDt)
  reaction.D <- Diagonal(x = reaction.stat_unc^2)

  reaction.S <- sysCompHandler$map(reaction.expDt, origSysDt, ret.mat = TRUE)
  reaction.P <- sysCompHandler$cov(updSysDt, ret.mat = TRUE)
  reaction.d <- reaction.expDt[,DATA-FIT_MEAN]
  chisquare(reaction.d,reaction.D,reaction.S,reaction.P)
}

getNDFreaction <- function(reaction_string) {
  n_points <- length(optExpDt[REAC==reaction_string]$IDX)
  n_pars <- length(optParamDt[ADJUSTABLE==TRUE]$IDX)
  #n_points - n_pars
  n_points
}

getChi2reaction_with_model_unc <- function(reaction_string) {

  reaction.expDt <- copy(optExpDt[REAC==reaction_string])
  reaction.model_cov <- model_covariance[reaction.expDt$IDX,reaction.expDt$IDX]

  reaction.expDt[, IDX := seq_len(.N)]
  reaction.stat_unc <- getDt_UNC(reaction.expDt)
  reaction.D <- Diagonal(x = reaction.stat_unc^2)
  reaction.S <- sysCompHandler$map(reaction.expDt, origSysDt, ret.mat = TRUE)
  reaction.P <- sysCompHandler$cov(updSysDt, ret.mat = TRUE)
  reaction.exp_cov <- reaction.S %*% reaction.P %*% t(reaction.S) + reaction.D

  reaction.tot_cov <- reaction.model_cov + reaction.exp_cov

  reaction.d <- reaction.expDt[,DATA-FIT_MEAN]

  # this should be the proper way? but cov is singular so it can't be solved
  as.vector(t(reaction.d)%*%solve(reaction.tot_cov,reaction.d))

}

chi2dt <- data.table(reaction=reactions,
  chisquare=lapply(reactions,getChi2reaction),
  chisquare_mod=lapply(reactions,getChi2reaction_with_model_unc),
  ndf=lapply(reactions,getNDFreaction))
# -------------------------------------------

# create model grid for experimental data
modSubents <- lapply(reactions, createSubentStub, en=energyGridrandomFiles)
modDt_post <- exforHandler$extractData(modSubents, ret.values=FALSE)
Smod <- exforHandler$getJac(modDt_post, RandomFilesNeedsDt, modSubents)
setkey(modDt_post, IDX, DIDX)

# for plotting uncertainty bands
tmpDt <- copy(RandomFilesNeedsDt)
setkey(tmpDt, IDX)
uncinfo <- cov.wt(t(sampledResults))
tmpDt[, V1:=uncinfo$center]
tmpDt[, UNC:=sqrt(diag(uncinfo$cov))]
postDt <- copy(modDt_post)

postDt[, MEAN:=as.vector(Smod %*% uncinfo$center)]
postDt[, UNC:=as.vector(sqrt(diag(Smod %*% diag(tmpDt$UNC^2) %*% t(Smod))))]
postDt[, MODE:=as.vector(Smod %*% allResults[,1])]
# --------------------------------------

for (curReac in reactions) {

    curExpDt <- optExpDt[REAC == curReac]
    curModDt <- modDt[REAC == curReac]
    curPostDt <- postDt[REAC == curReac]

    # label with chi2 information
    Chi2_exp_str <- paste0("$\\Chi^2_{exp}$ / n = ",format(chi2dt[reaction==curReac]$chisquare,digit=3)," / ",chi2dt[reaction==curReac]$ndf)
    Chi2_fit_str <- paste0("$\\Chi^2_{exp+fit}$ / n = ",format(chi2dt[reaction==curReac]$chisquare_mod,digit=3)," / ",chi2dt[reaction==curReac]$ndf)
    tex_label <- TeX(paste0(Chi2_exp_str," | ",Chi2_fit_str))
    
    ggp <- ggplot(curExpDt,aes(x = L1, y = DATA)) + theme_bw()
    ggp <- ggp + theme(axis.text=element_text(size=9),
                       axis.title=element_text(size=10),
                       plot.title=element_text(size=12),
                       plot.subtitle=element_text(size=7))
    ggp <- ggp + guides(col = "none")
    ggp <- ggp + xlab("energy (MeV)") + ylab("cross section (mbarn)")
    ggp <- ggp + labs(title=curReac,subtitle=tex_label)

    ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = "black",
                               linewidth = 0.2, width = 0.25)
    ggp <- ggp + geom_point(aes(x = L1, y = DATA), size=0.25)

    # plot the default TALYS model
    ggp <- ggp + geom_line(data=curModDt, aes(x = L1, y = DATA), col="black", linewidth=0.2, alpha=0.5,linetype = "dashed")
    #ggp <- ggp + geom_point(data=curModDt[,c("L1","DATA")], aes(x = L1, y = DATA), col="red", size=0.2)

    # plot the model posterior
    ggp <- ggp + geom_line(aes(x=L1, y=MEAN), data=curPostDt, col="green", linewidth=0.2)
    ggp <- ggp + geom_ribbon(aes(x=L1, ymin=MEAN-UNC, ymax=MEAN+UNC, y=MEAN), data=curPostDt,fill="green",alpha=0.3)
    ggp <- ggp + geom_line(aes(x=L1, y=MODE), data=curPostDt, col="red", linewidth=0.2)
    #ggp <- ggp + xlim(c(0,10))
    #ggp <- ggp + facet_wrap(~REAC, scales='free_y')


    #print(ggp)
    dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
    filepath <- file.path(plotPath, paste0('posterior_TALYS_with_gp_obs_', curReac,'.png'))
    ggsave(filepath, ggp, width = 8.65, height = 5.6, units = "cm", dpi = 300)
    filepath <- file.path(plotPath, paste0('posterior_TALYS_with_gp_obs_', curReac,'_LE.png'))
    ggsave(filepath, ggp+xlim(c(1,3)), width = 8.65, height = 5.6, units = "cm", dpi = 300)
}
