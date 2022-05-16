#################################################
#       SCRIPT Setup
##################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  source("./config/config.R")
  stop("No config file supplied, using default file config.R", call.=FALSE)
} else if (length(args) > 1) {
  stop("Script only accepts one argument.", call.=FALSE)
} else {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}

library(ggplot2)

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

origSysDt <- read_object(4, "origSysDt")
updSysDt <- read_object(4, "updSysDt")
expDt <- read_object(3, "expDt")

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


reactions <- expDt[,unique(REAC)]

for (curReac in reactions) {

    curExpDt <- expDt[REAC == curReac]
    ggp <- ggplot(curExpDt) + theme_bw()
    ggp <- ggp + theme(axis.text=element_text(size=9),
                       axis.title=element_text(size=10),
                       plot.title=element_text(size=12))
    ggp <- ggp + guides(col = "none")
    ggp <- ggp + xlab("energy [MeV]") + ylab("cross section [mbarn]")
    ggp <- ggp + ggtitle(curReac)

    ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = "green",
                               size = 0.5, width = 0.2)
    ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - ORIGUNC, ymax = DATA + ORIGUNC, col = EXPID),
                               size = 0.5, width = 0.3)
    ggp <- ggp + geom_point(aes(x = L1, y = DATA, col = EXPID))

    print(ggp)
    dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
    filepath <- file.path(plotPath, paste0('MLO_correction_', curReac,'.png'))
    ggsave(filepath, ggp, width = 8.65, height = 5.6, units = "cm", dpi = 300)
}

