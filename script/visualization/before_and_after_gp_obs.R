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

expDt <- read_object(3, "expDt", outdata_path=outdataPath)
optRes_before <- read_object(7, "optRes", outdata_path=outdataPath)
optRes_after <- read_object(10, "optRes", outdata_path=outdataPath)
updSysDt <- read_object(4, "updSysDt", outdata_path=outdataPath)
sysDt <- updSysDt[grepl('^EXPID-', EXPID)]
sysDt[, IDX:=seq_len(.N)]

# reconstruct experimental covariance matrix
sysCompHandler <- createSysCompHandler()
normHandler <- createSysCompNormHandler(dataref='DATAREF')
normHandler$addSysUnc('EXPID', unique(expDt$EXPID), 0, 0, TRUE)
sysCompHandler$addHandler(normHandler)
S_syserr <- sysCompHandler$map(expDt, sysDt, ret.mat=TRUE)
Cov_syserr <- sysCompHandler$cov(sysDt, ret.mat=TRUE)
totunc_exp <- sqrt(expDt$UNC^2 + rowSums((S_syserr %*% Cov_syserr) * S_syserr))

fitResult_before <- as.numeric(as.vector(optRes_before$fn))
fitUnc_before <- as.numeric(as.vector(optRes_before$stdAsyFitErrLM))

fitResult_after <- as.numeric(as.vector(optRes_after$fn))
fitUnc_after <- as.numeric(as.vector(optRes_after$stdAsyFitErrLM))

tmpExpDt <- copy(expDt)
setkey(tmpExpDt, IDX)
tmpExpDt[, UNC := totunc_exp]

tmpExpDt[, fit1 := fitResult_before]
tmpExpDt[, fitunc1 := fitUnc_before]

tmpExpDt[, fit2 := fitResult_after]
tmpExpDt[, fitunc2 := fitUnc_after]

ggp <- ggplot(data=tmpExpDt)
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=4),
                   axis.title=element_text(size=4),
                   strip.text=element_text(size=3))
ggp <- ggp + xlab('enegy [MeV]') + ylab('cross section [mbarn]')
# overlay experimental data
ggp <- ggp + geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC, col=EXPID), data=tmpExpDt)
ggp <- ggp + geom_point(aes(x=L1, y=DATA, col=EXPID), data=tmpExpDt, size=0.2)
# plot the models

ggp <- ggp + geom_ribbon(aes(x=L1, ymin=fit1-fitunc1, ymax=fit1+fitunc1), alpha=0.3, color="red",linewidth=0.2)
#ggp <- ggp + geom_errorbar(aes(x=L1, ymin=fit1-fitunc1, ymax=fit1+fitunc1), color="red", data=tmpExpDt,linewidth=0.1)


ggp <- ggp + geom_ribbon(aes(x=L1, ymin=fit2-fitunc2, ymax=fit2+fitunc2), alpha=0.3, color="green",linewidth=0.2)
#ggp <- ggp + geom_errorbar(aes(x=L1, ymin=fit2-fitunc2, ymax=fit2+fitunc2), color="green", data=tmpExpDt,linewidth=0.1)

ggp <- ggp + facet_wrap(~REAC, scales='free_y')

ggsave(file.path(plotPath, 'before_and_after_gp_obs.png'), ggp, units='cm', width=29.7/3, height=21/3)

ggp
