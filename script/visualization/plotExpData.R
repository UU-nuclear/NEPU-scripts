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

expDt <- read_object(3, "expDt")

reactions <- expDt[,unique(REAC)]

for (curReac in reactions) {

    curExpDt <- expDt[REAC == curReac]

    ggp <- ggplot(curExpDt,aes(x = L1, y = DATA)) + theme_bw()
    ggp <- ggp + scale_x_continuous(breaks=seq(0,30,5))
    ggp <- ggp + theme(axis.text=element_text(size=9),
                       axis.title=element_text(size=10),
                       plot.title=element_text(size=12),
                       plot.subtitle=element_text(size=7))
    ggp <- ggp + guides(col = "none")
    ggp <- ggp + xlab("energy (MeV)") + ylab("cross section (mbarn)")
    ggp <- ggp + labs(title=curReac)

    ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - UNC, ymax = DATA + UNC), col = curExpDt$EXPID,
                               size = 0.2, width = 0.25)
    ggp <- ggp + geom_point(aes(x = L1, y = DATA), size=0.25, col = curExpDt$EXPID)

    print(ggp)
    dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
    filepath <- file.path(plotPath, paste0('ExpData_', curReac,'.png'))
    ggsave(filepath, ggp, width = 8.65, height = 5.6, units = "cm", dpi = 300)
}