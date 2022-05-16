setwd("/home/joape785/eval-fe56")
source("config.R")
#source("config/config_2-30_fullParset_gp_obs.R")
#source("config/config_2-30_fullParset_gp_obs_2.R")
#source("config/config_2-55_fullParset_gp_obs.R")
#source("config/config_2-55_fullParset_gp_obs_reduced_tot.R")
#source("config/config_unmodified_2-55.R")
#scriptnr <- 8L
library(ggplot2)
library(stringr)


plot_path <- file.path(rootpath, outdataDir, "plots")


origSysDt <- read_object(4, "origSysDt")
updSysDt <- read_object(4, "updSysDt")
expDt <- read_object(3, "expDt")
refParamDt <- read_object(2, "refParamDt")

optRes <- read_object(7, "optRes")
P0 <- read_object(7, "P0")
X <- read_object(7, "X")
SX <- read_object(7, "S0")

fitResult <- as.numeric(as.vector(optRes$fn))
J <- optRes$jac
finalPars <- optRes$par
finalParCovmat <- optRes$parCovLM


fit_unc <- as.vector(optRes$stdAsyFitErrLM)

expDt[, FITUNC := fit_unc]
expDt[, XSECTFIT := fitResult]


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
expDt[, ORIGUNC := origUnc]
expDt[, UPDUNC := updUnc]

#D <- sparseMatrix(i=seq_along(statUnc),j=seq_along(statUnc), x=statUnc, dims=c(length(seq_along(statUnc)), length(seq_along(statUnc))))
D = Diagonal(x = statUnc^2)
P <- updX
updCov <- diag(statUnc^2) + S %*% updX %*% t(S)
d <- expDt[,DATA-XSECTFIT]
updCovInv <- inv_cov(D, S, updX)
#updCovInv <- solve(updCov)
redchi2tot <- t(d) %*% updCovInv %*% d / nrow(expDt)
redchi2tot <- chisquare(d, D, S, P, cholZ = NULL) / nrow(expDt)
redchi2tot_diag <- t(d) %*% diag(1/updUnc) %*% d / nrow(expDt)

# speeds up inversions of the experimental covariance matrix
cholZ <- makeCholZ(D, SX, X)

dpriorRef <- finalPars - 1
LpriorRef <- as.vector(crossprod(dpriorRef, solve(P0, dpriorRef)))
Lref1 <- as.vector(mult_xt_invCov_x(d, D, SX, X, cholZ = cholZ)) + LpriorRef
Lref2 <- as.vector(mult_xt_invCov_x(expDt[,DATA-DATAREF], D, SX, X, cholZ = cholZ)) + LpriorRef


residuals_stat <- (expDt$DATA - expDt$XSECTFIT) / expDt$UNC
residuals_sys_stat <- (expDt$DATA - expDt$XSECTFIT) / expDt$ORIGUNC
residuals_sys_upd_stat <- (expDt$DATA - expDt$XSECTFIT) / sqrt(expDt$UPDUNC^2 + expDt$FITUNC^2)

expDt[, NORMRESID := residuals_sys_upd_stat]
#hist( (expDt$DATA - expDt$XSECTFIT) / expDt$UNC, breaks=100)
#hist( (expDt$DATA - expDt$XSECTFIT) / expDt$ORIGUNC, breaks=100)


library(MASS)
library(actuar)
library(fitdistrplus)
library(moments)
#library(tikzDevice)

fit.norm <- fitdist(residuals_sys_upd_stat, "cauchy")  # we assume my_data ~ Normal(?,?)
result <- gofstat(f=fit.norm)
result$chisq
result$chisqdf
result$chisqpvalue
mean(residuals_sys_upd_stat)
sd(residuals_sys_upd_stat)
sum(residuals_sys_upd_stat*residuals_sys_upd_stat)/length(residuals_sys_upd_stat)

#tikz('normal.tex', standAlone = TRUE, width=5, height=5)
h <- hist( residuals_sys_upd_stat, pch=20, breaks=400, prob=TRUE, main="")
curve(dcauchy(x, fit.norm$estimate[1], fit.norm$estimate[2]), col="red", lwd=2, add=T)
#curve(pcauchy(x, fit.cauchy$estimate[1], fit.cauchy$estimate[2]), col="orange", lwd=2, add=T)
#dev.off()
Ntot <- length(expDt$NORMRESID)
avr <- mean(expDt$NORMRESID)
std <- sd(expDt$NORMRESID)
redChi2 <- sum(expDt$NORMRESID*expDt$NORMRESID)/(length(expDt$NORMRESID)-length(finalPars))
skw <- skewness(expDt$NORMRESID)
krt <- kurtosis(expDt$NORMRESID)
lableString <- sprintf("N: %1.0f\n Avr: %1.2f\n Std: %1.2f\n Skw: %1.2f\n Ex krt: %1.2f",Ntot, avr, std, skw, krt-3)

histPlot <- ggplot(expDt, aes(x=NORMRESID)) + theme_bw() + theme(aspect.ratio = 0.5) + 
  xlab(expression((y-f)/sigma)) + ylab("Density") + ggtitle("Normalized residuals") +
  geom_histogram(aes(y=..density..), color="black", fill="white") +
  annotate("text", x = Inf, y = Inf, label = "Preliminary", 
           hjust =1.2, vjust = 0.5, size=20, fontface="bold.italic", color="gray", alpha= 0.3, angle = 45) +
  #stat_function(fun = dnorm, args =  list(fit.norm$estimate[1], fit.norm$estimate[2]), aes(color="")) +
  annotate("text", x = Inf, y = Inf, label = lableString, hjust =1.5, vjust = 1.5, size=5)  #+
#labs(colour = "Fit")
print(histPlot)

## KS-plot normal
ggplot(expDt, aes(NORMRESID)) +
  stat_ecdf() + 
  stat_function(fun = pnorm, colour = "red")

ed <- ecdf(expDt$NORMRESID)
maxdiffidx <- which.max(abs(ed(expDt$NORMRESID)-pnorm(expDt$NORMRESID)))
maxdiffat <- expDt$NORMRESID[maxdiffidx]

ggplot(expDt, aes(NORMRESID)) +
  stat_ecdf() + 
  stat_function(fun = pnorm, colour = "red") + 
  geom_vline(xintercept=maxdiffat, lty=2)

## KS-plot cauchy
ggplot(expDt, aes(NORMRESID)) +
  stat_ecdf() + 
  stat_function(fun = pcauchy, colour = "red")

ed <- ecdf(expDt$NORMRESID)
maxdiffidx <- which.max(abs(ed(expDt$NORMRESID)-pcauchy(expDt$NORMRESID)))
maxdiffat <- expDt$NORMRESID[maxdiffidx]

ggplot(expDt, aes(NORMRESID)) +
  stat_ecdf() + 
  stat_function(fun = pcauchy, colour = "red") + 
  geom_vline(xintercept=maxdiffat, lty=2)


ggplot(expDt, aes(sample=NORMRESID)) + stat_qq() + stat_qq_line()



reactions <- unique(expDt$REAC)

for (curReac in reactions) {
  curExpDt <- expDt[REAC == curReac]
  avr <- mean(curExpDt$NORMRESID)
  std <- sd(curExpDt$NORMRESID)
  redChi2 <- sum(curExpDt$NORMRESID*curExpDt$NORMRESID)/length(curExpDt$NORMRESID)
  lableString <- sprintf("Mean: %1.2f \n Std: %1.2f \n Chi2/ndof:  %1.2f  ", avr,std, redChi2)
  
  
  ggp <- ggplot(curExpDt) + theme_bw()
  ggp <- ggp + xlab("Energy, Lab (MeV)") + ylab("Cross section (mbarn)")
  ggp <- ggp + ggtitle(curReac)
  
  ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = "black", size = 0.4)
  ggp <- ggp + geom_errorbar(aes(x = L1, ymin = DATA - ORIGUNC, ymax = DATA + ORIGUNC, col = EXPID), size = 0.2)
  ggp <- ggp + geom_point(aes(x = L1, y = DATA, col = EXPID), size = 0.5)
  #ggp <- ggp + geom_ribbon(aes(x = L1, ymin = XSECTFIT - 3*FITUNC, ymax = XSECTFIT + 3*FITUNC, alpha=0.5), fill = "grey70") 
  #ggp <- ggp + geom_ribbon(aes(x = L1, ymin = XSECTFIT - 2*FITUNC, ymax = XSECTFIT + 2*FITUNC, alpha=0.6), fill = "grey70") 
  ggp <- ggp + geom_ribbon(aes(x = L1, ymin = XSECTFIT - FITUNC, ymax = XSECTFIT + FITUNC, alpha=0.7), fill = "grey70") 
  ggp <- ggp + geom_line(aes(x = L1, y = XSECTFIT), col = "black",  size = 0.5)   
  ggp <- ggp + geom_line(aes(x = L1, y = DATAREF), col = "orange",  size = 0.5)   
  
  ggp <- ggp + scale_alpha(guide = 'none')
  ggp <- ggp + annotate("text", x = Inf, y = Inf, label = "Preliminary", 
                        hjust =1.2, vjust = 0.5, size=20, fontface="bold.italic", color="gray", alpha= 0.3, angle = 45)
  #ggp <- ggp + annotate("text", x = Inf, y = Inf, label = lableString, 
  #                     hjust =1.1, vjust = 1.1, size=5)  
  #+scale_color_manual(name = "LM-fit", values = c("fit" = "black"))
  # ggp <- ggp + scale_fill_continuous(name = "Data") 
  # ggp <- ggp
  
  print(ggp)
  #filepath <- file.path(plot_path, "gp_obs_fit") 
  #filename <- paste0(curReac, "_gp_obs_fit.pdf")
  #if(!dir.exists(filepath)) dir.create(file.path(filepath), showWarnings = FALSE, recursive = TRUE)
  #ggsave(file.path(filepath, filename), ggp, width = 36, height = 24, units = "cm", dpi = 300)  
  
}
