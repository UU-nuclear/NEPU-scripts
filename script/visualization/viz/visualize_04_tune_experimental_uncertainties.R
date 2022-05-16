source("config.R")
library(ggplot2)

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

library(ggplot2)



expDt[,{
    curExpDt <- .SD
    curReac <- .BY
    
    ggp <- ggplot(curExpDt) + theme_bw()
    ggp <- ggp + theme(axis.text=element_text(size=9),
                       axis.title=element_text(size=10),
                       plot.title=element_text(size=12))
    ggp <- ggp + guides(col = FALSE)
    ggp <- ggp + xlab("energy [MeV]") + ylab("cross section [mbarn]")
    ggp <- ggp + ggtitle(curReac)

    ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = EXPID,
                               size = 0.5, width = 0.2, alpha = 0.5)
    ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - ORIGUNC, ymax = DATA + ORIGUNC), col = "gray",
                               size = 0.6, width = 0.3)
    ggp <- ggp + geom_point(aes(x = L1, y = DATA), col = "black", size=0.5)
    ggp <- ggp + geom_line(aes(x = L1, y = DATAREF), size=1.5, color="white")
    ggp <- ggp + geom_line(aes(x = L1, y = DATAREF), size=1, color="black")
    print(ggp)
    dir.create(file.path(plotPath, "MLO_correction"), recursive=TRUE, showWarnings=FALSE)
    filename <- paste0("MLO_correction_",curReac, ".pdf")
    filepath <- file.path(plotPath, "MLO_correction", filename)
    ggsave(filepath, ggp, width = 40, height = 20, units = "cm", dpi = 300)
    }, by=REAC]

expDt[,{
  curExpDt <- .SD
  curReac <- .BY
  
  ggp <- ggplot(curExpDt) + theme_bw()
  ggp <- ggp + theme(axis.text=element_text(size=15),
                     axis.title=element_text(size=20),
                     plot.title=element_text(size=20))
  ggp <- ggp + guides(col = FALSE)
  ggp <- ggp + xlab("energy [MeV]") + ylab("cross section [mbarn]")
  ggp <- ggp + ggtitle(curReac)
  
  ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - ORIGUNC, ymax = DATA + ORIGUNC, col = EXPID),
                             size = 0.5, width = 0.2, alpha = 0.5)
  ggp <- ggp + geom_point(aes(x = L1, y = DATA, col = EXPID), size=0.5)
  #ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - ORIGUNC, ymax = DATA + ORIGUNC), col = "gray",
  #                           size = 0.6, width = 0.3)
  #ggp <- ggp + geom_point(aes(x = L1, y = DATA), col = "black", size=0.5)
  ggp <- ggp + guides(col = guide_legend(nrow = 20, byrow = TRUE))

  print(ggp)
  dir.create(file.path(plotPath, "data"), recursive=TRUE, showWarnings=FALSE)
  filename <- paste0("data_",curReac, ".pdf")
  filepath <- file.path(plotPath, "data", filename)
  ggsave(filepath, ggp, width = 40, height = 20, units = "cm", dpi = 300)
}, by=REAC]

expDt[,{
  curExpDt <- .SD
  curReac <- .BY
  
  ggp <- ggplot(curExpDt) + theme_bw()
  ggp <- ggp + theme(axis.text=element_text(size=15),
                     axis.title=element_text(size=20),
                     plot.title=element_text(size=20))
  ggp <- ggp + guides(col = FALSE)
  ggp <- ggp + xlab("energy [MeV]") + ylab("cross section [mbarn]")
  ggp <- ggp + ggtitle(curReac)
  
  ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - ORIGUNC, ymax = DATA + ORIGUNC, col = EXPID),
                             size = 0.5, width = 0.2, alpha = 0.5)
  ggp <- ggp + geom_point(aes(x = L1, y = DATA, col = EXPID), size=0.5)
  
  ggp <- ggp + geom_line(aes(x = L1, y = DATAREF), size=1)
  #ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - ORIGUNC, ymax = DATA + ORIGUNC), col = "gray",
  #                           size = 0.6, width = 0.3)
  #ggp <- ggp + geom_point(aes(x = L1, y = DATA), col = "black", size=0.5)
  ggp <- ggp + guides(col = guide_legend(nrow = 20, byrow = TRUE))
  
  print(ggp)
  dir.create(file.path(plotPath, "data_ref"), recursive=TRUE, showWarnings=FALSE)
  filename <- paste0("data_ref",curReac, ".pdf")
  filepath <- file.path(plotPath, "data_ref", filename)
  ggsave(filepath, ggp, width = 40, height = 20, units = "cm", dpi = 300)
}, by=REAC]
