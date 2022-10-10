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
tmpExpDt[, LMUNC := fitUnc]

D = Diagonal(x = getDt_UNC(expDt)^2)
chi2_ref <- chisquare(tmpExpDt[,DATA-DATAREF], D, SX, X)
chi2_lm <- chisquare(tmpExpDt[,DATA-LMFIT], D, SX, X)

ggp <- ggplot(data=tmpExpDt)
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=9),
                   axis.title=element_text(size=10),
                   strip.text=element_text(size=8))
ggp <- ggp + xlab('enegy [MeV]') + ylab('cross section [mbarn]')
# plot the model
ggp <- ggp + geom_line(aes(x=L1, y=LMFIT),col="red")
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=LMFIT-LMUNC, ymax=LMFIT+LMUNC),col="red",fill="red",alpha=0.5)
ggp <- ggp + facet_wrap(~REAC, scales='free_y')


#dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
ggsave(file.path(plotPath, 'posteriorXS.png'), ggp, units='cm', width=17.8, height=10)

