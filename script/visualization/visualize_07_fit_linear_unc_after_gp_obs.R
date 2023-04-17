#################################################
#       SCRIPT Setup
##################################################

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
extNeedsDt <- read_object(2, "extNeedsDt", outdata_path=outdataPathRun)
expDt <- read_object(3, "expDt", outdata_path=outdataPathRun)
updSysDt <- read_object(4, "updSysDt", outdata_path=outdataPathRun)
sysDt <- updSysDt[grepl('^EXPID-', EXPID)]
sysDt[, IDX:=seq_len(.N)]
needsDt <- read_object(1, "needsDt", outdata_path=outdataPathRun)
optRes <- read_object(10, "optRes", outdata_path=outdataPathRun)
P0 <- read_object(10, "P0", outdata_path=outdataPathRun)
X <- read_object(10, "X", outdata_path=outdataPathRun)
SX <- read_object(10, "S0", outdata_path=outdataPathRun)

fitResult <- as.numeric(as.vector(optRes$fn))
fitUnc <- as.numeric(as.vector(optRes$stdAsyFitErrLM))
J <- optRes$jac

############################
# reconstruct experimental covariance matrix
normHandler2 <- createSysCompNormHandler("DATAREF")
normHandler2$addSysUnc("EXPID", "", 0, 0, TRUE)

sysCompHandler2 <- createSysCompHandler()
sysCompHandler2$addHandler(normHandler2)

S <- sysCompHandler2$map(expDt, updSysDt, ret.mat = TRUE)
updX <- sysCompHandler2$cov(updSysDt, ret.mat = TRUE)
statUnc <- getDt_UNC(expDt)

totunc_exp <- sqrt(statUnc^2 + diag(S %*% updX %*% t(S)))
############################

tmpExpDt <- copy(expDt)
setkey(tmpExpDt, IDX)
tmpExpDt[, UNC := totunc_exp]
tmpExpDt[, LMFIT := fitResult]
tmpExpDt[, LMUNC := fitUnc]

D = Diagonal(x = getDt_UNC(expDt)^2)
chi2_ref <- chisquare(tmpExpDt[,DATA-DATAREF], D, SX, X)
chi2_lm <- chisquare(tmpExpDt[,DATA-LMFIT], D, SX, X)

# plot the results

#hist(tmpExpDt[,(DATA-DATAREF)/UNC], breaks=100)
#qqnorm(tmpExpDt[,DATA-DATAREF], pch = 1, frame = FALSE)
#qqline(tmpExpDt[,DATA-DATAREF], col = "steelblue", lwd = 2)
##avr(tmpExpDt[,DATA-DATAREF])
#all.moments(tmpExpDt[,(DATA-DATAREF)/UNC], order.max = 4, central = FALSE, absolute = FALSE, na.rm = FALSE)
#hist(tmpExpDt[,(DATA-LMFIT)/UNC], breaks=100)
#qqnorm(tmpExpDt[,DATA-LMFIT], pch = 1, frame = FALSE)
#qqline(tmpExpDt[,DATA-LMFIT], col = "steelblue", lwd = 2)
#all.moments(tmpExpDt[,(DATA-LMFIT)/UNC], order.max = 4, central = FALSE, absolute = FALSE, na.rm = FALSE)
#crossprod(tmpExpDt[,DATA-LMFIT], solve(D + SX %*% X %*% t(SX), tmpExpDt[,DATA-LMFIT]))

ggp <- ggplot(data=tmpExpDt)
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=4),
                   axis.title=element_text(size=4),
                   strip.text=element_text(size=3))
ggp <- ggp + xlab('energy (MeV)') + ylab('cross section (mbarn)')
# overlay experimental data
ggp <- ggp + geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC, col=EXPID), data=tmpExpDt, linewidth=0.2)
ggp <- ggp + geom_point(aes(x=L1, y=DATA, col=EXPID), data=tmpExpDt, size=0.1)
# plot the model
ggp <- ggp + geom_line(aes(x=L1, y=DATAREF), color="red",linewidth=0.2)
#ggp <- ggp + geom_point(aes(x=L1, y=DATAREF), color="red", size=0.2)
ggp <- ggp + geom_line(aes(x=L1, y=LMFIT),linewidth=0.2)
#ggp <- ggp + geom_point(aes(x=L1, y=LMFIT), size=0.2)
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=LMFIT-LMUNC, ymax=LMFIT+LMUNC), alpha=0.3)
#ggp <- ggp + geom_point(aes(x=L1, y=LMFIT), size=0.01)
#ggp <- ggp + geom_errorbar(aes(x=L1, ymin=LMFIT-LMUNC, ymax=LMFIT+LMUNC), alpha=0.2)
#ggp <- ggp + geom_ribbon(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC), alpha=0.3)
ggp <- ggp + facet_wrap(~REAC, scales='free')

dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
ggsave(file.path(plotPath, 'plot_posterior_xs_after_gp_obs.png'), ggp, units='cm', width=29.7/3, height=21/3)
