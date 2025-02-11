

library(hetGP)
library(MASS)
library(ggplot2)
library(data.table)

expDt <- read_object(3,"expDt")
modDt <- read_object(3, "modDt")

# test different GP covariance functions on the difficult data set
# 13840002
experiment <- "13840002"
curExpDt <- expDt[EXPID==experiment]
#curExpDt <- curExpDt[L1<3] # just to test the script on my laptop

reac <- "(24-CR-52(N,TOT),,SIG)"
curModDt <- modDt[REAC==reac]


model_hetGP_exp_SqrExp <- read_object(21,"model_hetGP_exp_SqrExp")
model_hetGP_exp_Mat3_2 <- read_object(21,"model_hetGP_exp_Mat3_2")
model_hetGP_exp_Mat5_2 <- read_object(21,"model_hetGP_exp_Mat5_2")

pred_hetGP_exp_SqrExp <- read_object(21,"pred_hetGP_exp_SqrExp")
pred_hetGP_exp_Mat3_2 <- read_object(21,"pred_hetGP_exp_Mat3_2")
pred_hetGP_exp_Mat5_2 <- read_object(21,"pred_hetGP_exp_Mat5_2")

pred_hetGP_exp_SqrExp[,KERNEL:="Sqr. Exp."]
pred_hetGP_exp_Mat3_2[,KERNEL:="Matern 3/2"]
pred_hetGP_exp_Mat5_2[,KERNEL:="Matern 5/2"]

GP_predictions <- rbind(pred_hetGP_exp_SqrExp,pred_hetGP_exp_Mat3_2,pred_hetGP_exp_Mat5_2)

ggp <- ggplot(data=curExpDt)
ggp <- ggp + theme_bw() #+ theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=4),
                   axis.title=element_text(size=4),
                   strip.text=element_text(size=3))
ggp <- ggp + geom_line(aes(x=L1, y=DATA))
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=mean-sqrt(nugs), ymax=mean+sqrt(nugs), col=KERNEL),data=GP_predictions, alpha=0., linetype = "dashed")
ggp <- ggp + geom_line(aes(x=L1, y=mean, col=KERNEL),data=GP_predictions)
print(ggp)

# It appears that I cannot fit the talys model with the Matern3/2 kernel, it is probably not well sutied for that
# Consequently I get strange results when trying to fit the data with the estimated length scale (although the length-scale
# appears to be reasonable: 2.86 MeV for the main process and 19.9 for the noise process). The Matern3/2 kernel fits the
# resonance data instead and its mean does not appear smooth, like the Matern5/2. Thereby the estimated nuggets does not
# model the variance of data around the smooth trend, as intended.

# The Matern5/2 and the Square Exponential kernel results are very close to eachother

compareGP(model_hetGP_exp_SqrExp,model_hetGP_exp_Mat5_2)


ggp <- ggplot(data=curExpDt)
ggp <- ggp + theme_bw() #+ theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=10),
                   axis.title=element_text(size=15),
                   plot.title=element_text(size=12),
                   plot.subtitle=element_text(size=7))
ggp <- ggp + xlab("energy (MeV)") + ylab("cross section (mbarn)")
ggp <- ggp + geom_line(aes(x=L1, y=DATA))
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=mean-sqrt(nugs), ymax=mean+sqrt(nugs)),data=GP_predictions[KERNEL=="Matern 5/2"], alpha=0.2, linetype = "dashed", col="green", fill="green")
ggp <- ggp + geom_line(aes(x=L1, y=mean),data=GP_predictions[KERNEL=="Matern 5/2"], col="green")
ggp <- ggp + geom_line(aes(x=L1, y=DATA),data=curModDt, col="red")
ggp <- ggp + xlim(1,5)
print(ggp)


plotPath <- "/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/NEPU-scripts/pics-for-JEFF-meeting-2023-spring"

filepath <- file.path(plotPath, paste0('hetGPfit.png'))
ggsave(filepath, ggp, width = 3*8.65, height = 3*5.6, units = "cm", dpi = 300)
