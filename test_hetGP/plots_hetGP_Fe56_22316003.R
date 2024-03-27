# FOr the figure from WONDER-2023: source("config/config-Fe56.R")


args = commandArgs(trailingOnly=TRUE)

if(length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
} else if (length(args) > 1) {
  stop("Script only accepts one argument.", call.=FALSE)
}

scriptnr <- 21L
overwrite <- FALSE

dir.create(file.path(outdataPath,scriptnr), recursive=TRUE, showWarnings=FALSE)

library(hetGP)
library(MASS)
library(ggplot2)
library(data.table)

expDt <- read_object(3,"expDt")
modDt <- read_object(3, "modDt")

experiment <- "22316003"
curExpDt <- expDt[EXPID==experiment]
#curExpDt <- curExpDt[L1<3] # just to test the script on my laptop

reac <- curExpDt[,unique(REAC)]
curModDt <- modDt[REAC==reac]

# scale to barn
curExpDt[,DATA:=DATA*1E-03]
curModDt[,DATA:=DATA*1E-03]

# pred_hetGP_talys_Mat5_2 <- read_object(21,"pred_hetGP_talys_Mat5_2")
# pred_hetGP_exp_Mat5_2 <- read_object(21,"pred_hetGP_exp_Mat5_2")
# pred_hetGP_talys_Mat5_2 <- read_object(21,"pred_hetGP_talys_SqrExp")
# pred_hetGP_exp_SqrExp <- read_object(21,"pred_hetGP_exp_SqrExp")

model_hetGP_exp_Mat5_2 <- read_object(21,"model_hetGP_exp_Mat5_2")

Xmin <- curExpDt[,min(L1)] - 0.1
Xmax <- 20#curExpDt[,max(L1)] + 0.3

Xgrid <- matrix(seq(Xmin,Xmax,by=0.1), ncol=1)
pred_hetGP <- as.data.table(predict(x = Xgrid, object = model_hetGP_exp_Mat5_2))
pred_hetGP[,L1:=Xgrid]

pred_hetGP[,mean:=mean*1E-03]
pred_hetGP[,nugs:=nugs*1E-06]


ggp <- ggplot(data=curExpDt)
ggp <- ggp + theme_bw() #+ theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=8),
                   axis.title=element_text(size=8),
                   strip.text=element_text(size=8))
#ggp <- ggp + geom_line(aes(x=L1, y=DATA))
ggp <- ggp + xlab("energy (MeV)") + ylab("cross section (barn)")
ggp <- ggp + geom_point(aes(x=L1, y=DATA),size=0.1)
ggp <- ggp + geom_line(aes(x=L1, y=DATA),data=curModDt[L1<21],col='red',linewidth=0.4)
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=mean-sqrt(nugs), ymax=mean+sqrt(nugs)),
                          data=pred_hetGP, alpha=0.125, linetype = "dashed",
                          fill='green', col='green', linewidth=0.35)
ggp <- ggp + geom_line(aes(x=L1, y=mean),data=pred_hetGP,col='green',linewidth=0.4)
ggp <- ggp + coord_cartesian(xlim=c(0,20))
print(ggp)

#ggp + scale_color_manual(breaks=c('GP', 'default','prior'),
#                       values=c('extra sys. unc.'='green', 'only stat. unc.'='red', 'orig. sys. unc.'='black','prior'='blue'))

plotPath <- "/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/pics-for-JEFF-meeting-2023-spring"

filepath <- file.path(plotPath, paste0('hetGPfit_Fe56_22316003.pdf'))
ggsave(filepath, ggp, width = 130*2/3, height = 56.1, units = "mm")

# ==========================================================
Xgrid <- matrix(curExpDt[,L1], ncol=1)
pred_hetGP <- as.data.table(predict(x = Xgrid, object = model_hetGP_exp_Mat5_2))
pred_hetGP[,L1:=Xgrid]

curExpDt[,GP_mean:=1.e-03*pred_hetGP[,mean]]
curExpDt[,GP_nugs:=1.e-06*pred_hetGP[,nugs]]

curExpDt[,RESIDUAL:=(DATA-GP_mean)/sqrt(GP_nugs)]
#curExpDt[,RESIDUAL:=(DATA-GP_mean)]


ggp <- ggplot(data=curExpDt,aes(x=RESIDUAL)) +
theme_bw() +
theme(axis.text=element_text(size=8),
                   axis.title=element_text(size=8),
                   strip.text=element_text(size=8)) +
xlab(TeX("normalized residual")) + 
geom_histogram(aes(y=after_stat(density)),binwidth=0.1,color="grey", fill="grey",linewidth=0.0) +
#geom_freqpoly(aes(y=after_stat(density)),binwidth=0.1,color="black", fill="grey",linewidth=0.01) +
#stat_bin(aes(y=after_stat(density)),binwidth=0.1,color="black", fill="grey",linewidth=0.01) +
coord_cartesian(xlim=c(-4,4)) +
stat_function(fun = dnorm, args = list(mean = 0, sd = 1),color='green',linetype = "dashed", linewidth=0.5,xlim=c(-4,4)) +
scale_y_continuous(position = "right") +
scale_x_continuous(breaks=seq(-4,4,by=2),minor_breaks=seq(-5,5))
ggp

filepath <- file.path(plotPath, paste0('hetGPfit_Fe56_22316003_residual.pdf'))
ggsave(filepath, ggp, width = 130/3, height = 56.1, units = "mm")

# ===========================================================

curExpDt[,ORIG_UNC:=1.e-03*ORIG_UNC]
# Make the same type of figure but for the reported uncertainties

ggp <- ggplot(data=curExpDt)
ggp <- ggp + theme_bw() #+ theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=8),
                   axis.title=element_text(size=8),
                   strip.text=element_text(size=8))
#ggp <- ggp + geom_line(aes(x=L1, y=DATA))
ggp <- ggp + xlab("energy (MeV)") + ylab("cross section (barn)")
#ggp <- ggp + geom_point(aes(x=L1, y=DATA),size=0.1)
ggp <- ggp + geom_errorbar(aes(x=L1, ymin=DATA-ORIG_UNC,ymax=DATA+ORIG_UNC),size=0.1)
ggp <- ggp + geom_line(aes(x=L1, y=DATA),data=curModDt[L1<21],col='red',linewidth=0.4)
ggp <- ggp + coord_cartesian(xlim=c(0,20))
print(ggp)

plotPath <- "/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/pics-for-JEFF-meeting-2023-spring"

filepath <- file.path(plotPath, paste0('hetGPfit_Fe56_22316003_noGP.pdf'))
ggsave(filepath, ggp, width = 130*2/3, height = 56.1, units = "mm")

curExpDt[,DATAREF:=1.e-03*DATAREF]
curExpDt[,RESIDUAL_TALYS:=(DATA-DATAREF)/sqrt(ORIG_UNC)]

ggp <- ggplot(data=curExpDt,aes(x=RESIDUAL_TALYS)) +
theme_bw() +
theme(axis.text=element_text(size=8),
                   axis.title=element_text(size=8),
                   strip.text=element_text(size=8)) +
xlab(TeX("normalized residual")) + 
geom_histogram(aes(y=after_stat(density)),binwidth=0.1,color="grey", fill="grey",linewidth=0.0) +
#geom_freqpoly(aes(y=after_stat(density)),binwidth=0.1,color="black", fill="grey",linewidth=0.01) +
#stat_bin(aes(y=after_stat(density)),binwidth=0.1,color="black", fill="grey",linewidth=0.01) +
coord_cartesian(xlim=c(-10,10)) +
stat_function(fun = dnorm, args = list(mean = 0, sd = 1),color='red',linetype = "dashed", linewidth=0.5,xlim=c(-10,10)) +
scale_y_continuous(position = "right") +
scale_x_continuous(breaks=seq(-10,10,by=5),minor_breaks=seq(-10,10,by=2.5))
ggp

filepath <- file.path(plotPath, paste0('hetGPfit_Fe56_22316003_residual_noGP.pdf'))
ggsave(filepath, ggp, width = 130/3, height = 56.1, units = "mm")