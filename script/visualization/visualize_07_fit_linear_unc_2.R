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
optRes <- read_object(7, "optRes", outdata_path=outdataPathRun)
P0 <- read_object(7, "P0", outdata_path=outdataPathRun)
X <- read_object(7, "X", outdata_path=outdataPathRun)
SX <- read_object(7, "S0", outdata_path=outdataPathRun)

fitResult <- as.numeric(as.vector(optRes$fn))
J <- optRes$jac


# plot the results

# reconstruct experimental covariance matrix
sysCompHandler <- createSysCompHandler()
normHandler <- createSysCompNormHandler(dataref='DATAREF')
normHandler$addSysUnc('EXPID', unique(expDt$EXPID), 0, 0, TRUE)
sysCompHandler$addHandler(normHandler)
S_syserr <- sysCompHandler$map(expDt, sysDt, ret.mat=TRUE)
Cov_syserr <- sysCompHandler$cov(sysDt, ret.mat=TRUE)
totunc_exp <- sqrt(expDt$UNC^2 + rowSums((S_syserr %*% Cov_syserr) * S_syserr))
tmpExpDt <- copy(expDt)
setkey(tmpExpDt, IDX)
tmpExpDt[, UNC := totunc_exp]
tmpExpDt[, LMFIT := fitResult]

D = Diagonal(x = getDt_UNC(expDt)^2)
chi2_ref <- chisquare(tmpExpDt[,DATA-DATAREF], D, SX, X)
chi2_lm <- chisquare(tmpExpDt[,DATA-LMFIT], D, SX, X)

hist(tmpExpDt[,(DATA-DATAREF)/UNC], breaks=100)
qqnorm(tmpExpDt[,DATA-DATAREF], pch = 1, frame = FALSE)
qqline(tmpExpDt[,DATA-DATAREF], col = "steelblue", lwd = 2)
avr(tmpExpDt[,DATA-DATAREF])
all.moments(tmpExpDt[,(DATA-DATAREF)/UNC], order.max = 4, central = FALSE, absolute = FALSE, na.rm = FALSE)
hist(tmpExpDt[,(DATA-LMFIT)/UNC], breaks=100)
qqnorm(tmpExpDt[,DATA-LMFIT], pch = 1, frame = FALSE)
qqline(tmpExpDt[,DATA-LMFIT], col = "steelblue", lwd = 2)
all.moments(tmpExpDt[,(DATA-LMFIT)/UNC], order.max = 4, central = FALSE, absolute = FALSE, na.rm = FALSE)
#crossprod(tmpExpDt[,DATA-LMFIT], solve(D + SX %*% X %*% t(SX), tmpExpDt[,DATA-LMFIT]))

plot_lm <- function(){
  ggp <- ggplot(data=tmpExpDt)
  ggp <- ggp + theme_bw() + theme(legend.position="none")
  ggp <- ggp + theme(axis.text=element_text(size=9),
                     axis.title=element_text(size=10),
                     strip.text=element_text(size=8))
  ggp <- ggp + xlab('enegy [MeV]') + ylab('cross section [mbarn]')
  # overlay experimental data
  ggp <- ggp + geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC, col=EXPID), data=tmpExpDt)
  ggp <- ggp + geom_point(aes(x=L1, y=DATA, col=EXPID), data=tmpExpDt, size=0.2)
  # plot the model
  ggp <- ggp + geom_line(aes(x=L1, y=DATAREF), color="red")
    ggp <- ggp + geom_line(aes(x=L1, y=LMFIT))
  #ggp <- ggp + geom_ribbon(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC), alpha=0.3)
  ggp <- ggp + facet_wrap(~REAC, scales='free_y')
  
  ggp  
}

print(plot_lm)

#dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
#ggsave(file.path(plotPath, 'plot_posterior_xs.png'), ggp,
#       units='cm', width=17.8, height=10)



