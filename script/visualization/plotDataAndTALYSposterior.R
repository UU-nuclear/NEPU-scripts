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

extNeedsDt <- read_object(2, "extNeedsDt")
origSysDt <- read_object(4, "origSysDt")
updSysDt <- read_object(4, "updSysDt")
expDt <- read_object(3, "expDt")
modDt <- read_object(3,"modDt")
allResults <- read_object(9, 'allResults')

optParamDt <- read_object(7,"optParamDt")
optRes <- read_object(7, "optRes")
P0 <- read_object(7, "P0")
X <- read_object(7, "X")
SX <- read_object(7, "S0")

reactions <- expDt[,unique(REAC)]

normHandler <- createSysCompNormHandler("DATAREF")
normHandler$addSysUnc("EXPID", "", 0, 0, TRUE)

sysCompHandler <- createSysCompHandler()
sysCompHandler$addHandler(normHandler)

S <- sysCompHandler$map(expDt, origSysDt, ret.mat = TRUE)
origX <- sysCompHandler$cov(origSysDt, ret.mat = TRUE)
updX <- sysCompHandler$cov(updSysDt, ret.mat = TRUE)
statUnc <- getDt_UNC(expDt)

origUnc <- sqrt(statUnc^2 + diag(S %*% origX %*% t(S))) 
updUnc <- sqrt(statUnc^2 + diag(S %*% updX %*% t(S)))

setkey(expDt, IDX)
expDt[, ORIGUNC := origUnc]
expDt[, UPDUNC := updUnc]

# ----------Calculate the Chi-squares -----------
fitResult <- as.numeric(as.vector(optRes$fn))
J <- optRes$jac
finalPars <- optRes$par
finalParCovmat <- optRes$parCovLM
fit_unc <- as.vector(optRes$stdAsyFitErrLM)

expDt[, FITUNC := fit_unc]
expDt[, XSECTFIT := fitResult]

D <- Diagonal(x = statUnc^2)
P <- updX
d <- expDt[,DATA-XSECTFIT]
chi2_tot <- chisquare(d,D,S,P)

getChi2reaction <- function(reaction_string) {
  normHandler <- createSysCompNormHandler("DATAREF")
  normHandler$addSysUnc("EXPID", "", 0, 0, TRUE)

  sysCompHandler <- createSysCompHandler()
  sysCompHandler$addHandler(normHandler)

  reaction.expDt <- expDt[REAC==reaction_string]
  reaction.expDt[, IDX := seq_len(.N)]
  reaction.stat_unc <- getDt_UNC(reaction.expDt)
  reaction.D <- Diagonal(x = reaction.stat_unc^2)

  reaction.S <- sysCompHandler$map(reaction.expDt, origSysDt, ret.mat = TRUE)
  reaction.P <- sysCompHandler$cov(updSysDt, ret.mat = TRUE)
  reaction.d <- reaction.expDt[,DATA-XSECTFIT]
  chisquare(reaction.d,reaction.D,reaction.S,reaction.P)
}

getNDFreaction <- function(reaction_string) {
  n_points <- length(expDt[REAC==reaction_string]$IDX)
  n_pars <- length(optParamDt[ADJUSTABLE==TRUE]$IDX)
  #n_points - n_pars
  n_points
}

chi2dt <- data.table(reaction=reactions,chisquare=lapply(reactions,getChi2reaction),ndf=lapply(reactions,getNDFreaction))
# -------------------------------------------

# create model grid for experimental data
modSubents <- lapply(reactions, createSubentStub, en=energyGrid)
#modSubents <- lapply(reactions, createSubentStub, en=energyGridrandomFiles)
modDt_post <- exforHandler$extractData(modSubents, ret.values=FALSE)
Smod <- exforHandler$getJac(modDt_post, extNeedsDt, modSubents)
setkey(modDt_post, IDX, DIDX)

# for plotting uncertainty bands
tmpDt <- copy(extNeedsDt)
setkey(tmpDt, IDX)
uncinfo <- cov.wt(t(allResults))
tmpDt[, V1:=uncinfo$center]
tmpDt[, UNC:=sqrt(diag(uncinfo$cov))]
postDt <- copy(modDt_post)

postDt[, DATA:=as.vector(Smod %*% tmpDt$V1)]
postDt[, UNC:=as.vector(diag(Smod %*% diag(tmpDt$UNC) %*% t(Smod)))]
# --------------------------------------

for (curReac in reactions) {

    curExpDt <- expDt[REAC == curReac]
    curModDt <- modDt[REAC == curReac]
    curPostDt <- postDt[REAC == curReac]

    # label with chi2 information
    tex_label <- TeX(paste0("$\\Chi^2$ / n = ",format(chi2dt[reaction==curReac]$chisquare,digit=3)," / ",chi2dt[reaction==curReac]$ndf))

    ggp <- ggplot(curExpDt,aes(x = L1, y = DATA)) + theme_bw()
    ggp <- ggp + scale_x_continuous(breaks=seq(0,30,5))
    ggp <- ggp + theme(axis.text=element_text(size=9),
                       axis.title=element_text(size=10),
                       plot.title=element_text(size=12),
                       plot.subtitle=element_text(size=7))
    ggp <- ggp + guides(col = "none")
    ggp <- ggp + xlab("energy (MeV)") + ylab("cross section (mbarn)")
    ggp <- ggp + labs(title=curReac,subtitle=tex_label)

    ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = "black",
                               size = 0.2, width = 0.25)
    ggp <- ggp + geom_point(aes(x = L1, y = DATA), size=0.25)

    # plot the default TALYS model
    ggp <- ggp + geom_line(data=curModDt[,c("L1","DATA")], aes(x = L1, y = DATA), col="red", size=0.2)
    ggp <- ggp + geom_point(data=curModDt[,c("L1","DATA")], aes(x = L1, y = DATA), col="red", size=0.2)

    # plot the model posterior
    ggp <- ggp + geom_line(aes(x=L1, y=DATA), data=curPostDt, col="green", size=0.2)
    ggp <- ggp + geom_ribbon(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC), alpha=0.3, data=curPostDt,fill="green")
    #ggp <- ggp + facet_wrap(~REAC, scales='free_y')


    print(ggp)
    dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
    filepath <- file.path(plotPath, paste0('posterior_TALYS_', curReac,'.png'))
    ggsave(filepath, ggp, width = 8.65, height = 5.6, units = "cm", dpi = 300)
}