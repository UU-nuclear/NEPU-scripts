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
optRes <- read_object(7, "optRes", outdata_path=outdataPathRun)
P0 <- read_object(7, "P0", outdata_path=outdataPathRun)
X <- read_object(7, "X", outdata_path=outdataPathRun)
SX <- read_object(7, "S0", outdata_path=outdataPathRun)

#optRes <- read_object(10, "optRes", outdata_path=outdataPathRun)
#P0 <- read_object(10, "P0", outdata_path=outdataPathRun)
#X <- read_object(10, "X", outdata_path=outdataPathRun)
#SX <- read_object(10, "S0", outdata_path=outdataPathRun)

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

tmpExpDt[, RESIDUAL := (DATA - LMFIT)]

ggp <- ggplot(data=tmpExpDt)
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=9),
                   axis.title=element_text(size=10),
                   strip.text=element_text(size=8))
ggp <- ggp + xlab('enegy [MeV]') + ylab('residual [mbarn]')

ggp <- ggp + geom_errorbar(aes(x=L1, ymin=RESIDUAL-UNC, ymax=RESIDUAL+UNC, col=EXPID), data=tmpExpDt, alpha=0.3)
ggp <- ggp + geom_point(aes(x=L1, y=RESIDUAL), data=tmpExpDt, size=0.2, alpha=0.6)

ggp <- ggp + facet_wrap(~REAC, scales='free_y')
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=-LMUNC, ymax=+LMUNC, col='red'), alpha=0.3)

# make histograms of the residuals

# (n,tot)-channel
totExpDt <- tmpExpDt[REAC=="(26-FE-56(N,TOT),,SIG)"]

histplot <- ggplot(data=totExpDt) + geom_histogram(aes(x=(DATA - LMFIT), color=EXPID, fill=EXPID),binwidth=50, position="identity",alpha=0.3)

histplot <- ggplot(data=totExpDt)
histplot <- histplot + geom_histogram(data = totExpDt[EXPID==22316003],aes(x=(DATA - LMFIT), color=EXPID, fill=EXPID),binwidth=50, position="identity",alpha=0.3)
histplot <- histplot + geom_histogram(data = totExpDt[EXPID==13764002],aes(x=(DATA - LMFIT), color=EXPID, fill=EXPID),binwidth=50, position="identity",alpha=0.3)
histplot <- histplot + geom_histogram(data = totExpDt[EXPID==41325003],aes(x=(DATA - LMFIT), color=EXPID, fill=EXPID),binwidth=50, position="identity",alpha=0.3)
histplot <- histplot + geom_histogram(data = totExpDt[EXPID==10037005],aes(x=(DATA - LMFIT), color=EXPID, fill=EXPID),binwidth=50, position="identity",alpha=0.3)

histplot2 <- ggplot(data=tmpExpDt) + geom_histogram(aes(x=(DATA - LMFIT), color=REAC, fill=REAC),binwidth=50, position="identity",alpha=0.3)

ggplot(data=tmpExpDt) + geom_histogram(aes(x=log(((DATA - LMFIT)/UNC)**2),y=..density.., color=REAC, fill=REAC),binwidth=0.5, position="identity",alpha=0.3)

#ggplot(data=tmpExpDt[REAC==reactions[3]]) + geom_histogram(aes(x=log((DATA - LMFIT)**2/(UNC**2+122.34**2))),binwidth=0.5, position="identity",alpha=0.3)

chi2_vals <- seq(-10,5,by=0.1)
chi2_density <- dchisq(exp(chi2_vals),1)

tmp_dt <- data.table(chi2_vals)
tmp_dt[,density := chi2_density]

# plot the GP observable

# read in the results of 07_5_addGPobs.R

optExpDt <- read_object(7,"optExpDt")
optSysDt <- read_object(7,"optSysDt")
optGpDt <- read_object(7,"optGpDt")
mapAssignment <- read_object(7,"mapAssignment")
reacHandlerGPobs <- read_object(7,"reacHandlerGPobs")
sysCompHandler <- read_object(7,"sysCompHandler")
Xk <- read_object(7,"Xk")
S0k <- read_object(7,"S0k")


#sysErrDt <- optSysDt[grepl("EXPID-",optSysDt$EXPID)]
sysErrDt <- optSysDt[optSysDt[!grepl("TALYS",EXPID)]]
#sysErrDt <- optSysDt
sysErrDt[, IDX := seq_len(.N)]
setkey(sysErrDt, IDX)
S <- sysCompHandler$map(optExpDt,sysErrDt, ret.mat = TRUE)
U <- sysCompHandler$cov(sysErrDt, ret.mat = TRUE)

sys_pwl_Dt <- optSysDt[ERRTYPE=="pw"]
sys_pwl_Dt[, IDX := seq_len(.N)]
setkey(sys_pwl_Dt, IDX)

S_pwl <- sysCompHandler$map(optExpDt, sys_pwl_Dt, ret.mat = TRUE)
#P_pwl <- sysCompHandler$cov(sys_pwl_Dt, ret.mat = TRUE)
P_pwl <- sysCompHandler$cov(sys_pwl_Dt, optGpDt, ret.mat = TRUE)


GPmean <- t(S_pwl %*% P_pwl) %*% mult_invCov_x(optExpDt$DATA-optExpDt$REFDATA,Diagonal(x=optExpDt$UNC^2),S,U)
sys_pwl_Dt[,DATA:=as.vector(GPmean)]

currReac <- "(26-FE-56(N,P)25-MN-56,,SIG)"
#currReac <- "(26-FE-56(N,TOT),,SIG)"


currReac2 <- mapAssignment[REAC==currReac]$EXPID

ggp <- ggplot(data=tmpExpDt[REAC==currReac])
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=9),
                   axis.title=element_text(size=10),
                   strip.text=element_text(size=8))
ggp <- ggp + xlab('enegy [MeV]') + ylab('residual [mbarn]')

ggp <- ggp + geom_errorbar(aes(x=L1, ymin=RESIDUAL-UNC, ymax=RESIDUAL+UNC, col=EXPID), alpha=0.3)
ggp <- ggp + geom_point(aes(x=L1, y=RESIDUAL), size=0.2, alpha=0.6)
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=-LMUNC, ymax=+LMUNC, col='red'), alpha=0.3)
ggp <- ggp + geom_line(data=sys_pwl_Dt[EXPID==currReac2],aes(x=EN,y=DATA))

# Try it for a single reaction =========================

optExpDt[, TOTUNC := totunc_exp]
optExpDt[, LMFIT := fitResult]
optExpDt[, LMUNC := fitUnc]
optExpDt[, RESIDUAL := DATA-LMFIT]

#curREAC <- "(26-FE-56(N,2N)26-FE-55,,SIG)"
curREAC <- "(26-FE-56(N,P)25-MN-56,,SIG)"



"(26-FE-56(N,P)25-MN-56,,SIG)"
"(26-FE-56(N,EL)26-FE-56,,SIG)" 
"(26-FE-56(N,TOT),,SIG)"
"(26-FE-56(N,INL)26-FE-56,,SIG)"
"(26-FE-56(N,D)25-MN-55,,SIG)"
"(26-FE-56(N,A)24-CR-53,,SIG)"  
"(26-FE-56(N,2N)26-FE-55,,SIG)"
"(26-FE-56(N,T)25-MN-54,,SIG)"  
"(26-FE-56(N,N+P)25-MN-55,,SIG)"

curExpIds <- expDt[REAC == curREAC, unique(EXPID)]

curExpDt <- optExpDt[EXPID %in% curExpIds,]
curExpDt[, IDX := seq_len(.N)]

# create sysDt only containing systematic errors of
# experiments in current reaction channel
curSysDt <- optSysDt[EXPID %in% paste0("EXPID-",curExpIds) | 
                  EXPID==mapAssignment[REAC==curREAC]$EXPID,]   
curSysDt[, IDX := seq_len(.N)]


currGpDt <- optGpDt[EXPID==mapAssignment[REAC==curREAC]$EXPID]

sys_pwl_Dt <- curSysDt[ERRTYPE=="pw"]
sys_pwl_Dt[, IDX := seq_len(.N)]
setkey(sys_pwl_Dt, IDX)


S <- sysCompHandler$map(curExpDt,curSysDt, ret.mat = TRUE)
#U <- sysCompHandler$cov(curSysDt, ret.mat = TRUE)
U <- sysCompHandler$cov(curSysDt, currGpDt, ret.mat = TRUE)

S_pwl <- sysCompHandler$map(curExpDt, sys_pwl_Dt, ret.mat = TRUE)
P_pwl <- sysCompHandler$cov(sys_pwl_Dt, currGpDt, ret.mat = TRUE)


Sigma22inv_x_residual <- mult_invCov_x(curExpDt$DATA-curExpDt$REFDATA,Diagonal(x=curExpDt$UNC^2),S,U)

# extract the mean of the gp
# given by K(x,X) (K(X,X)+SIGMA_EXP)^-1 residual
GPmean <- t(S_pwl %*% P_pwl) %*% Sigma22inv_x_residual

# extract the uncertainty of the gp
# the covariance is K(x,x) -  K(x,X) (K(X,X)+SIGMA_EXP)^-1 K(X,x)
# the diagonal of K(x,x) is just the sigma hyper-parameter squared
GPunc <-  rep((currGpDt[PARNAME=="sigma"]$PARVAL)**2,nrow(sys_pwl_Dt)) - 
  diag(mult_xt_invCov_x(S_pwl %*% P_pwl,Diagonal(x=curExpDt$UNC^2),S,U))
GPunc <- sqrt(GPunc)

sys_pwl_Dt[,DATA:=as.vector(GPmean)]
sys_pwl_Dt[,UNC:=as.vector(GPunc)]

ggp <- ggplot(data=curExpDt)
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=9),
                   axis.title=element_text(size=10),
                   strip.text=element_text(size=8))
ggp <- ggp + xlab('enegy [MeV]') + ylab('residual [mbarn]')

ggp <- ggp + geom_errorbar(aes(x=L1, ymin=RESIDUAL-TOTUNC, ymax=RESIDUAL+TOTUNC, col=EXPID), alpha=0.3)
ggp <- ggp + geom_point(aes(x=L1, y=RESIDUAL), size=0.2, alpha=0.6)
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=-LMUNC, ymax=+LMUNC), col='red', fill='red', alpha=0.3)

ggp <- ggp + geom_line(data=sys_pwl_Dt,aes(x=EN,y=DATA))
ggp <- ggp + geom_ribbon(data=sys_pwl_Dt,aes(x=EN, ymin=DATA-UNC, ymax=DATA+UNC), alpha=0.3)