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
#extNeedsDt <- read_object(2, "extNeedsDt")
origSysDt <- read_object(4, "origSysDt")
updSysDt <- read_object(4, "updSysDt")
expDt <- read_object(3, "expDt")
modDt <- read_object(3,"modDt")
optExpDt <- read_object(7, "optExpDt")
allResults <- read_object(12, 'allResults')

optParamDt <- read_object(10,"optParamDt")
optRes <- read_object(10, "optRes")
P0 <- read_object(10, "P0")
X <- read_object(10, "X")
SX <- read_object(10, "S0")

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
modSubents <- lapply(reactions, createSubentStub, en=energyGridrandomFiles)
modDt_post <- exforHandler$extractData(modSubents, ret.values=FALSE)
setkey(modDt_post, IDX, DIDX)

uncinfo <- cov.wt(t(allResults))

modDt_post[,V1:=uncinfo$center]
modDt_post[, UNC:=sqrt(diag(uncinfo$cov))]

ggp <- ggplot(expDt) + theme_bw()
ggp <- ggp + scale_x_continuous(breaks=seq(0,200,10))
ggp <- ggp + theme(axis.text=element_text(size=9),
                   axis.title=element_text(size=10),
                   plot.title=element_text(size=12),
                   plot.subtitle=element_text(size=7))
ggp <- ggp + guides(col = "none")
ggp <- ggp + xlab("energy (MeV)") + ylab("cross section (mbarn)")

ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = "black",
                           size = 0.2, width = 0.25)
ggp <- ggp + geom_point(aes(x = L1, y = DATA), size=0.25)

# plot the model posterior
ggp <- ggp + geom_line(data=modDt_post, aes(x=L1, y=V1), col="green", size=0.2)
ggp <- ggp + geom_ribbon(data=modDt_post, aes(x=L1, ymin=V1-UNC, ymax=V1+UNC), alpha=0.3, fill="green")
ggp <- ggp + facet_wrap(~REAC, scales='free_y')

#for (curReac in reactions) {
#    curExpDt <- expDt[REAC == curReac]
#    curPostDt <- modDt_post[REAC == curReac]
#
#    ggp <- ggplot(curExpDt,aes(x = L1, y = DATA)) + theme_bw()
#    ggp <- ggp + scale_x_continuous(breaks=seq(0,30,5))
#    ggp <- ggp + theme(axis.text=element_text(size=9),
#                       axis.title=element_text(size=10),
#                       plot.title=element_text(size=12),
#                       plot.subtitle=element_text(size=7))
#    ggp <- ggp + guides(col = "none")
#    ggp <- ggp + xlab("energy (MeV)") + ylab("cross section (mbarn)")
#    ggp <- ggp + labs(title=curReac,subtitle=tex_label)
#
#    ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = "black",
#                               size = 0.2, width = 0.25)
#    ggp <- ggp + geom_point(aes(x = L1, y = DATA), size=0.25)
#
#    # plot the model posterior
#    ggp <- ggp + geom_line(aes(x=L1, y=DATA), data=curPostDt, col="green", size=0.2)
#    ggp <- ggp + geom_ribbon(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC), alpha=0.3, data=curPostDt,fill="green")
#}