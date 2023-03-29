#################################################
#       SCRIPT Setup
##################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}

library(ggplot2)
#library(latex2exp)

expDt <- read_object(1, "expDtFull")

reactions <- expDt[,unique(REAC)]

for (curReac in reactions) {

    curExpDt <- expDt[REAC == curReac]

    ggp <- ggplot(curExpDt,aes(x = L1, y = DATA, col=EXPID)) + theme_bw()
    #ggp <- ggp + scale_x_continuous(breaks=seq(0,30,5))
    ggp <- ggp + theme(axis.text=element_text(size=9),
                       axis.title=element_text(size=10),
                       plot.title=element_text(size=12),
                       plot.subtitle=element_text(size=7),
                       legend.text = element_text(size=5))
    #ggp <- ggp + guides(col = "none")
    ggp <- ggp + xlab("energy (MeV)") + ylab("cross section (mbarn)")
    #ggp <- ggp + labs(title=curReac)

    #ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - UNC, ymax = DATA + UNC), col = curExpDt$EXPID,
    #                           linewidth = 0.2, width = 0.25)
    ggp <- ggp + geom_point(aes(x = L1, y = DATA, col=EXPID), size=0.25)

    ggp <- ggp + scale_x_continuous(breaks=energyGridrandomFiles, limits=c(0.0,30))
#
#    energy_cutoff <- 1.8
#    ## Fit  a piece-wise linear model on the talys-energy grid
#    lm1 <- lm(formula = DATA ~ bs(L1, df = NULL, knots = energyGridrandomFiles[energyGridrandomFiles>=energy_cutoff & energyGridrandomFiles<25], degree = 1),
#              data    = curExpDt[L1>energy_cutoff])
#
#    lm2 <- lm(formula = DATA ~ bs(L1, df = NULL, knots = energyGridrandomFiles[energyGridrandomFiles<25], degree = 1),
#              data    = curExpDt)
#
#    lm3 <- lm(formula = DATA ~ bs(L1, df = NULL, knots = energyGridrandomFiles[energyGridrandomFiles>=energy_cutoff & energyGridrandomFiles<25], degree = 1),
#              data    = curExpDt)
#
#    ## Create a data frames to hold prediction
#    newdat1 <- data.frame(L1 = energyGridrandomFiles[energyGridrandomFiles>=energy_cutoff])
#    newdat1$DATA <- predict(lm1, newdata = newdat1)
#
#    newdat2 <- data.frame(L1 = energyGridrandomFiles)
#    newdat2$DATA <- predict(lm2, newdata = newdat2)
#
#    newdat3 <- data.frame(L1 = energyGridrandomFiles[energyGridrandomFiles>=energy_cutoff])
#    newdat3$DATA <- predict(lm3, newdata = newdat3)
#
#    ## Plot the previous plot with a regression line
#    ggp <- ggp +
#      geom_line(data=newdat1,aes(x = L1, y = DATA),inherit.aes=FALSE,col='red') +
#      geom_point(data=newdat1,aes(x = L1, y = DATA),inherit.aes=FALSE,col='red') +
#      geom_line(data=newdat2,aes(x = L1, y = DATA),inherit.aes=FALSE,col='green') +
#      geom_point(data=newdat2,aes(x = L1, y = DATA),inherit.aes=FALSE,col='green') +
#      geom_line(data=newdat3,aes(x = L1, y = DATA),inherit.aes=FALSE,col='blue') +
#      geom_point(data=newdat3,aes(x = L1, y = DATA),inherit.aes=FALSE,col='blue') +
#      ylim(c(0,10000))
#
#    print(ggp)

    dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
    filepath <- file.path(plotPath, paste0('ExpData_nounc_Full', curReac,'.png'))
    ggsave(filepath, ggp, width = 16, height = 9, units = "cm", dpi = 300)
}



