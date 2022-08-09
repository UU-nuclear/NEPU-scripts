#################################################
#       SCRIPT Setup
##################################################
args = commandArgs(trailingOnly=TRUE)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}

library(ggplot2)
library(moments)

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################


outdataPathRun <- outdataPath

expDt <- read_object(3, "expDt", outdata_path=outdataPathRun)
optRes_10 <- read_object(10, "optRes", outdata_path=outdataPathRun)
optRes_7 <- read_object(7, "optRes", outdata_path=outdataPathRun)

fitResult_10 <- as.numeric(as.vector(optRes$fn))
fitUnc_10 <- as.numeric(as.vector(optRes$stdAsyFitErrLM))
fitResult_7 <- as.numeric(as.vector(optRes$fn))
fitUnc_7 <- as.numeric(as.vector(optRes$stdAsyFitErrLM))

# plot the results
tmpExpDt <- copy(expDt)

setkey(tmpExpDt, IDX)
tmpExpDt[, UNC := totunc_exp]
tmpExpDt[, LMFIT_10 := fitResult_10]
tmpExpDt[, LMUNC_10 := fitUnc_10]
tmpExpDt[, LMFIT_7 := fitResult_7]
tmpExpDt[, LMUNC_7 := fitUnc_7]


ggp <- ggplot(data=tmpExpDt)
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=9),
                   axis.title=element_text(size=10),
                   strip.text=element_text(size=8))
ggp <- ggp + xlab('enegy [MeV]') + ylab('cross section [mbarn]')
# overlay experimental data
ggp <- ggp + geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC, col=EXPID), data=tmpExpDt)
ggp <- ggp + geom_point(aes(x=L1, y=DATA, col=EXPID), data=tmpExpDt, size=0.2)

# plot the models
ggp <- ggp + geom_line(aes(x=L1, y=LMFIT_10),color="red")
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=LMFIT_10-LMUNC_10, ymax=LMFIT_10+LMUNC_10), alpha=0.3,color="red")
ggp <- ggp + geom_errorbar(aes(x=L1, ymin=LMFIT_10-LMUNC_10, ymax=LMFIT_10+LMUNC_10), data=tmpExpDt,color="red")

ggp <- ggp + geom_line(aes(x=L1, y=LMFIT_7))
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=LMFIT_7-LMUNC_7, ymax=LMFIT_7+LMUNC_7), alpha=0.3)
ggp <- ggp + geom_errorbar(aes(x=L1, ymin=LMFIT_7-LMUNC_7, ymax=LMFIT_7+LMUNC_7), data=tmpExpDt)

ggp <- ggp + facet_wrap(~REAC, scales='free_y')

ggp  



#dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
#ggsave(file.path(plotPath, 'plot_posterior_xs_after_gp_obs.png'), ggp, units='cm', width=17.8, height=10)



