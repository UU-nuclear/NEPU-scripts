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
extNeedsDt <- read_object(2, "extNeedsDt")
origSysDt <- read_object(4, "origSysDt")
updSysDt <- read_object(4, "updSysDt")
expDt <- read_object(3, "expDt")
modDt <- read_object(3,"modDt")
allResults <- read_object(9, 'allResults')

optExpDt <- read_object(6, "optExpDt")
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

# create model grid for experimental data
modSubents <- lapply(reactions, createSubentStub, en=energyGrid)
#modSubents <- lapply(reactions, createSubentStub, en=energyGridrandomFiles)
modDt_post <- exforHandler$extractData(modSubents, ret.values=FALSE)
Smod <- exforHandler$getJac(modDt_post, extNeedsDt, modSubents)
setkey(modDt_post, IDX, DIDX)

# for plotting uncertainty bands
extNeedsDtMod <- extNeedsDt[,{
  stopifnot(all(L2 == 0) & all(L3 == 0))
  list(L1 = defineEnergyGrid(L1, energyGridrandomFiles, enPolicy="compgrid"),
       L2 = 0, L3 = 0, V1 = 0)
}, by=c("PROJECTILE", "ELEMENT", "MASS", "REAC")]

extNeedsDtMod[, IDX := seq_len(.N)]

sampledResults <- allResults[,2:ncol(allResults)]
optResult <- allResults[,1]
postDt <- copy(extNeedsDtMod)
setkey(tmpDt, IDX)
uncinfo <- cov.wt(t(sampledResults))
postDt[, OPT:=optResult]
postDt[, V1:=uncinfo$center]
postDt[, UNC:=sqrt(diag(uncinfo$cov))]
# --------------------------------------
postDt[REAC=="CS/TOT",REAC:="(26-FE-56(N,TOT),,SIG)"]
postDt[REAC=="CS/EL",REAC:="(26-FE-56(N,EL)26-FE-56,,SIG)"]
postDt[REAC=="CS/REAC/110000/TOT",REAC:="(26-FE-56(N,N+P)25-MN-55,,SIG)"]
postDt[REAC=="CS/REAC/001000/TOT",REAC:="(26-FE-56(N,D)25-MN-55,,SIG)"]
postDt[REAC=="CS/REAC/100000/TOT",REAC:="(26-FE-56(N,INL)26-FE-56,,SIG)"]
postDt[REAC=="CS/REAC/010000/TOT",REAC:="(26-FE-56(N,P)25-MN-56,,SIG)"]
postDt[REAC=="CS/REAC/200000/TOT",REAC:="(26-FE-56(N,2N)26-FE-55,,SIG)"]
postDt[REAC=="CS/REAC/000001/TOT",REAC:="(26-FE-56(N,A)24-CR-53,,SIG)"]
postDt[REAC=="CS/REAC/000100/TOT",REAC:="(26-FE-56(N,T)25-MN-54,,SIG)"]
# --------------------------------------



curExpDt <- expDt
curPostDt <- postDt

# label with chi2 information
ggp <- ggplot(curExpDt) + theme_bw()
ggp <- ggp + scale_x_continuous(breaks=seq(0,200,20))
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

# plot the model posterior
ggp <- ggp + geom_line(aes(x=L1, y=OPT), data=curPostDt, col="red", size=0.2)
ggp <- ggp + geom_line(aes(x=L1, y=V1), data=curPostDt, col="green", size=0.2)
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=V1-UNC, ymax=V1+UNC), alpha=0.3, data=curPostDt,fill="green")
ggp + facet_wrap(~REAC, scales='free_y')

# and the model MAP
#ggp <- ggp + geom_line(aes(x=L1, y=XSECTFIT_MAP), col="blue", size=0.2)


#print(ggp)
dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
filepath <- file.path(plotPath, paste0('posterior_TALYS_', curReac,'.png'))
ggsave(filepath, ggp, width = 8.65, height = 5.6, units = "cm", dpi = 300)
