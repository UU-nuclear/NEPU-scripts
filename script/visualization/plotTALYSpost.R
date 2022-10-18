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
library(ggnewscale)

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

subents <- read_object(1, "subents")
extNeedsDt <- read_object(2, "extNeedsDt")
origSysDt <- read_object(4, "origSysDt")
updSysDt <- read_object(4, "updSysDt")
expDt <- read_object(3, "expDt")
modDt <- read_object(3,"modDt")
#allResults <- read_object(9, 'allResults')
allResults <- read_object(12, 'allResults')

optExpDt <- read_object(6, "optExpDt")
# optParamDt <- read_object(7,"optParamDt")
# optRes <- read_object(7, "optRes")
# P0 <- read_object(7, "P0")
# X <- read_object(7, "X")
# SX <- read_object(7, "S0")
optParamDt <- read_object(10,"optParamDt")
optRes <- read_object(10, "optRes")
P0 <- read_object(10, "P0")
X <- read_object(10, "X")
SX <- read_object(10, "S0")

# found a possible bug in the parsing of the talys files
# when the talys cross section result is zero NA appears in the
# parsed output. For now I will just replace NAs with 0, and print a warning

if(any(is.na(allResults))) {
  cat("Warning: NAs appear in the sampled cross sections. \n
    A bug in the parsing of the talys files causes NAs to appear when TALYS predict 0. \n
    This bug should be fixed, but for now replacing NAs with 0.\n")
  allResults[is.na(allResults)] <- 0
}

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
setkey(postDt, IDX)
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

ggp <- ggplot(expDt) + theme_bw() +
  scale_x_continuous(breaks=seq(0,200,20)) +
  theme(axis.text=element_text(size=8),
                     axis.title=element_text(size=10),
                     plot.title=element_text(size=8),
                     plot.subtitle=element_text(size=7)) +
  guides(col = "none") +
  xlab("energy (MeV)") + ylab("cross section (mbarn)") +
  geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = "black",
                             size = 0.2, width = 0.25) +
  geom_point(aes(x = L1, y = DATA), size=0.25)

# plot the model posterior
ggp <- ggp + new_scale_colour() +
  geom_line(aes(x=L1, y=OPT, col='mode'), data=postDt, size=0.2) +
  geom_line(aes(x=L1, y=V1, col='mean'), data=postDt, size=0.2) +
  geom_ribbon(aes(x=L1, ymin=V1-UNC, ymax=V1+UNC), fill='green', alpha=0.3, data=postDt) +
  scale_color_manual(name='posterior',
                       breaks=c('mode', 'mean'),
                       values=c('mode'='red', 'mean'='green')) +
  facet_wrap(~REAC, scales='free_y')


#print(ggp)
dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
filepath <- file.path(plotPath, 'posterior_TALYS.png')
ggsave(filepath, ggp, width = 2*16, height = 2*9, units = "cm", dpi = 300)

# individual plots per channel
curPlotPath <- file.path(plotPath,'talysPosteriorPlots')
dir.create(curPlotPath, recursive=TRUE, showWarnings=FALSE)

for(curReac in unique(expDt$REAC))
{
  cur_expDt <- expDt[REAC==curReac]
  cur_postDt <- postDt[REAC==curReac]

  plot <- ggplot(cur_expDt) + theme_bw() +
    labs(title=curReac) +
    scale_x_continuous(breaks=seq(0,200,20)) +
    theme(axis.text=element_text(size=8),
                       axis.title=element_text(size=10),
                       plot.title=element_text(size=8),
                       plot.subtitle=element_text(size=7)) +
    guides(col = "none") +
    xlab("energy (MeV)") + ylab("cross section (mbarn)") +
    geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = "black",
                               size = 0.2, width = 0.25) +
    geom_point(aes(x = L1, y = DATA), size=0.25)

    # plot the model posterior
  plot <- plot + new_scale_colour() +
    geom_line(aes(x=L1, y=OPT, col='mode'), data=cur_postDt, size=0.2) +
    geom_line(aes(x=L1, y=V1, col='mean'), data=cur_postDt, size=0.2) +
    geom_ribbon(aes(x=L1, ymin=V1-UNC, ymax=V1+UNC), fill='green', alpha=0.3, data=cur_postDt) +
    scale_color_manual(name='posterior',
                         breaks=c('mode', 'mean'),
                         values=c('mode'='red', 'mean'='green'))

  curPlotPath <- file.path(plotPath,'talysPosteriorPlots',paste0('talys_',curReac,'.png'))
  ggsave(curPlotPath, plot, width = 16, height = 9, units = "cm", dpi = 300)
}

curReac <- "(26-FE-56(N,TOT),,SIG)"
plot <- ggplot(cur_expDt) + theme_bw() +
    labs(title=curReac) +
    theme(axis.text=element_text(size=8),
                       axis.title=element_text(size=10),
                       plot.title=element_text(size=8),
                       plot.subtitle=element_text(size=7)) +
    guides(col = "none") +
    xlab("energy (MeV)") + ylab("cross section (mbarn)") +
    #geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = "black",
    #                           size = 0.2, width = 0.25) +
    geom_point(aes(x = L1, y = DATA), size=0.25)

    # plot the model posterior
  plot <- plot + new_scale_colour() +
    geom_line(aes(x=L1, y=DATAREF, col='prior'), data=cur_expDt, size=0.2) +
    geom_line(aes(x=L1, y=OPT, col='mode'), data=cur_postDt, size=0.2) +
    geom_line(aes(x=L1, y=V1, col='mean'), data=cur_postDt, size=0.2) +
    geom_ribbon(aes(x=L1, ymin=V1-UNC, ymax=V1+UNC), fill='green', alpha=0.3, data=cur_postDt) +
    scale_color_manual(name='posterior',
                         breaks=c('mode', 'mean','prior'),
                         values=c('mode'='red', 'mean'='green', 'prior'='blue'))

  plot + xlim(1,3)
