#################################################
#       SCRIPT Setup
##################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args)==1) {
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
modDt <- read_object(3,"modDt")

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
    curModDt <- modDt[REAC == curReac]

    ggp <- ggplot(curExpDt,aes(x = L1, y = DATA)) + theme_bw()
    ggp <- ggp + scale_x_continuous(breaks=seq(0,30,2))
    ggp <- ggp + theme(axis.text=element_text(size=9),
                       axis.title=element_text(size=10),
                       plot.title=element_text(size=12))
    ggp <- ggp + guides(col = "none")
    ggp <- ggp + xlab("energy (MeV)") + ylab("cross section (mbarn)")
    ggp <- ggp + ggtitle(curReac)

    ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = "black",
                               size = 0.2, width = 0.25)
    ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - UNC, ymax = DATA + UNC), col = curExpDt$EXPID,
                               size = 0.2, width = 0.25)
    #ggp <- ggp + geom_point(aes(x = L1, y = DATA), size=0.25, col = curExpDt$EXPID)

    ggp <- ggp + geom_line(data=curModDt[,c("L1","DATA")],aes(x = L1, y = DATA),col = "black", size=0.2)

    print(ggp)
    dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
    filepath <- file.path(plotPath, paste0('reference_TALYS_', curReac,'.png'))
    ggsave(filepath, ggp, width = 8.65, height = 5.6, units = "cm", dpi = 300)
}